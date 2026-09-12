"""The index of a data root: what it reads, what it remembers, what it re-reads.

The point of the index is that a second visit costs nothing, so the test that
matters is the one showing that an unchanged file is not read again, and that a
file which was added IS.
"""
import numpy as np
from astropy.io import fits

from pca2d import scan


def write_tfits(path, instrument="NIRPS", snr=(100.0, 120.0, 140.0),
                exptime=300.0, mjd=60000.5):
    """A file with the headers a scan reads, and nothing else in it."""
    primary = fits.PrimaryHDU()
    primary.header["INSTRUME"] = instrument
    primary.header["EXPTIME"] = exptime
    name = "FluxAB" if instrument == "SPIROU" else "FluxA"
    science = fits.ImageHDU(data=np.zeros((2, 4), dtype=np.float32), name=name)
    for i, value in enumerate(snr):
        science.header["EXTSN%03d" % i] = value
    science.header["MJDMID"] = mjd
    fits.HDUList([primary, science]).writeto(path, overwrite=True)


def root_with(tmp_path, **objects):
    tmp_path.mkdir(parents=True, exist_ok=True)
    for name, n in objects.items():
        folder = tmp_path / name
        folder.mkdir()
        for i in range(n):
            write_tfits(str(folder / ("%04dt.fits" % i)),
                        instrument="SPIROU" if name == "TOI1452" else "NIRPS",
                        snr=(10.0 * (i + 1), 10.0 * (i + 2)),
                        exptime=60.0 * (i + 1), mjd=60000.0 + i)
    return str(tmp_path)


def test_one_file_gives_the_instrument_the_snr_and_the_exposure(tmp_path):
    path = str(tmp_path / "0001t.fits")
    write_tfits(path, snr=(100.0, 120.0, 140.0), exptime=240.0)
    record = scan.scan_file(path)
    assert record["instrument"] == "NIRPS"
    assert record["snr"] == 120.0, "the median of the per-order EXTSN"
    assert record["exptime"] == 240.0
    assert record["mjd"] == 60000.5


def test_the_index_lives_outside_the_data_root(tmp_path):
    """A data root can be read-only or shared, and is never written to."""
    root = root_with(tmp_path / "data", PROXIMA=2)
    home = str(tmp_path / "home")
    index, _tally = scan.update(root, index=scan.load(root, home), home=home)
    assert scan.save(index, root, home).startswith(home)
    assert sorted(p.name for p in (tmp_path / "data").rglob("*")) == \
        ["0000t.fits", "0001t.fits", "PROXIMA"], "nothing added under the root"
    # two roots whose last folder has the same name do not share one index
    other = scan.index_path(str(tmp_path / "elsewhere" / "data"), home)
    assert other != scan.index_path(root, home)


def test_a_second_visit_reads_only_what_was_added(tmp_path):
    root = root_with(tmp_path, PROXIMA=3, GJ1=2)
    home = str(tmp_path / "home")
    index, tally = scan.update(root, home=home)
    assert (tally["read"], tally["kept"]) == (5, 0)
    scan.save(index, root, home)

    index, tally = scan.update(root, index=scan.load(root, home), home=home)
    assert (tally["read"], tally["kept"]) == (0, 5), "nothing changed, nothing read"

    write_tfits(str(tmp_path / "GJ1" / "0009t.fits"), exptime=900.0)
    index, tally = scan.update(root, index=index, home=home)
    assert (tally["read"], tally["kept"]) == (1, 5)
    assert scan.summary(index, "GJ1")["files"] == 3

    (tmp_path / "GJ1" / "0009t.fits").unlink()
    index, tally = scan.update(root, index=index, home=home)
    assert tally["gone"] == 1 and scan.summary(index, "GJ1")["files"] == 2


def test_a_replaced_file_is_read_again(tmp_path):
    """Same name, different content: the size and the modification time say so."""
    root = root_with(tmp_path, PROXIMA=1)
    home = str(tmp_path / "home")
    index, _tally = scan.update(root, home=home)
    assert scan.summary(index, "PROXIMA")["snr"] == 15.0
    path = str(tmp_path / "PROXIMA" / "0000t.fits")
    write_tfits(path, snr=(500.0, 700.0), exptime=60.0)
    import os
    info = os.stat(path)
    os.utime(path, (info.st_atime, info.st_mtime + 10))
    index, tally = scan.update(root, index=index, home=home)
    assert tally["read"] == 1
    assert scan.summary(index, "PROXIMA")["snr"] == 600.0


def test_the_summary_is_what_a_row_shows(tmp_path):
    root = root_with(tmp_path, PROXIMA=3, TOI1452=2)
    index, _tally = scan.update(root, home=str(tmp_path / "home"))
    rows = {row["object"]: row for row in scan.summaries(index)}
    assert rows["PROXIMA"]["files"] == 3
    assert rows["PROXIMA"]["instrument"] == "NIRPS"
    assert rows["TOI1452"]["instrument"] == "SPIROU", \
        "read from the file, never from a configuration"
    assert rows["PROXIMA"]["snr"] == np.median([15.0, 25.0, 35.0])
    assert rows["PROXIMA"]["exptime"] == 120.0
    assert (rows["PROXIMA"]["first"], rows["PROXIMA"]["last"]) == (60000.0, 60002.0)
    assert scan.summary(index, "NOBODY")["files"] == 0


def test_an_unreadable_file_is_recorded_once_and_not_read_again(tmp_path):
    root = root_with(tmp_path, PROXIMA=1)
    (tmp_path / "PROXIMA" / "0002t.fits").write_text("not a FITS file")
    home = str(tmp_path / "home")
    index, tally = scan.update(root, home=home)
    assert tally["read"] == 2
    row = scan.summary(index, "PROXIMA")
    assert row["files"] == 2 and row["snr"] == 15.0, \
        "the unreadable one counts as a file and contributes no number"
    _index, tally = scan.update(root, index=index, home=home)
    assert tally["read"] == 0, "a file that cannot be read is not re-read at"\
        " every visit"


def test_an_index_from_another_version_is_rebuilt_rather_than_believed(tmp_path):
    root = root_with(tmp_path, PROXIMA=1)
    home = str(tmp_path / "home")
    index, _tally = scan.update(root, home=home)
    index["version"] = scan.VERSION + 1
    scan.save(index, root, home)
    assert scan.load(root, home)["objects"] == {}
