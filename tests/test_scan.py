"""The index of a data root: what it reads, what it remembers, what it re-reads.

The point of the index is that a second visit costs nothing, so the test that
matters is the one showing that an unchanged file is not read again, and that a
file which was added IS.
"""
import numpy as np
from astropy.io import fits

from pca2d import scan


def write_tfits(path, instrument="NIRPS", snr=(100.0, 120.0, 140.0),
                exptime=300.0, mjd=60000.5, mag=None):
    """A file with the headers a scan reads, and nothing else in it.

    The brightness goes in under the keyword that instrument writes, and those
    are not the same band: NIRPS gives J, SPIRou's OBJMAG is H.
    """
    primary = fits.PrimaryHDU()
    primary.header["INSTRUME"] = instrument
    primary.header["EXPTIME"] = exptime
    if mag is not None:
        if instrument == "SPIROU":
            primary.header["OBJMAG"] = mag
        else:
            primary.header["ESO OCS TARG JMAG"] = mag
            primary.header["ESO OCS TARG IMAG"] = 0.0   # a field not filled in
            primary.header["GUIMAGN"] = -9999.9         # a guider with no value
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


def test_the_names_come_before_the_numbers_and_each_object_is_announced(tmp_path):
    """A first scan of a campaign is minutes: the list cannot stay empty for it."""
    root = root_with(tmp_path, PROXIMA=3, GJ1=2)
    listed, done, files = {}, [], []
    scan.update(root, home=str(tmp_path / "home"), on_listed=listed.update,
                on_object=lambda name, k, n: done.append((name, k, n)),
                on_file=lambda name, i, n: files.append((name, i, n)))
    assert listed == {"GJ1": 2, "PROXIMA": 3}, \
        "the names and the counts cost one folder listing, before any header"
    assert done == [("GJ1", 2, 2), ("PROXIMA", 3, 3)], \
        "each object is announced when it is finished, known == total"
    assert files[0] == ("GJ1", 1, 2) and files[-1] == ("PROXIMA", 3, 3)


def test_the_numbers_arrive_every_ten_spectra_marked_as_estimates(tmp_path):
    """Ten spectra already give a campaign's SNR and exposure time: they go up
    at once, and the caller shows them as estimates until the rest are read."""
    root = root_with(tmp_path, PROXIMA=25)
    done = []
    index, _tally = scan.update(root, home=str(tmp_path / "home"), every=10,
                                on_object=lambda n, k, t: done.append((k, t)))
    assert done == [(10, 25), (20, 25), (25, 25)], \
        "every ten, then once complete"
    partial = [(k, t) for k, t in done if k < t]
    assert partial and all(k < t for k, t in partial), \
        "an estimate is recognisable by known < total"
    assert scan.summary(index, "PROXIMA")["files"] == 25


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


def test_the_first_files_read_are_spread_over_the_campaign(tmp_path):
    """Names sort by date, so the first ten files are the first night and their
    median is that night's weather. Measured on GJ 1: 134 against 163."""
    order = scan.spread(list(range(100)), first=10)
    assert sorted(order) == list(range(100)), "the same files, once each"
    assert order[:10] == [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
    assert len(set(order[:10])) == 10
    for prefix in (10, 20, 50):
        half = sum(1 for i in order[:prefix] if i >= 50)
        assert abs(half - prefix / 2) <= 1, \
            "any prefix covers both halves of the campaign"
    assert scan.spread([1, 2, 3], first=10) == [1, 2, 3], "nothing to spread"
    assert scan.spread([], first=10) == []

    # and a real scan reads in that order
    root = root_with(tmp_path, PROXIMA=30)
    seen = []
    scan.update(root, home=str(tmp_path / "home"), every=10,
                on_file=lambda name, i, n: seen.append(i))
    assert seen == list(range(1, 31))


def test_a_partial_median_is_of_the_campaign_not_of_its_first_night(tmp_path):
    """The SNR rises through the campaign; ten spread files must not report the
    beginning of it."""
    root = str(tmp_path)
    folder = tmp_path / "PROXIMA"
    folder.mkdir(parents=True)
    for i in range(100):
        write_tfits(str(folder / ("%04dt.fits" % i)), snr=(float(i + 1),),
                    exptime=60.0, mjd=60000.0 + i)
    seen = {}

    def on_object(name, known, total):
        if known == 10:
            seen["ten"] = scan.summary(index[0], name)["snr"]

    index = [scan.load(root, str(tmp_path / "home"))]
    index[0], _tally = scan.update(root, index=index[0], every=10,
                                   home=str(tmp_path / "home"),
                                   on_object=on_object)
    assert seen["ten"] == 46.0, "the median of 1, 11, 21, ... 91"
    assert scan.summary(index[0], "PROXIMA")["snr"] == 50.5, "all of them"
    # read in order, the first ten would have given 5.5, an order of magnitude
    # away from the campaign's value
    assert abs(seen["ten"] - 50.5) < 6.0


def test_the_magnitude_comes_with_the_band_it_is_in(tmp_path):
    """NIRPS writes J, SPIRou's OBJMAG is H (checked against SIMBAD on both
    SPIRou campaigns). One unlabelled column would be wrong by a magnitude."""
    nirps = str(tmp_path / "n.fits")
    spirou = str(tmp_path / "s.fits")
    write_tfits(nirps, instrument="NIRPS", mag=5.328)
    write_tfits(spirou, instrument="SPIROU", mag=10.452)
    assert scan.scan_file(nirps)["mag"] == 5.328
    assert scan.scan_file(nirps)["mag_band"] == "J"
    assert scan.scan_file(spirou)["mag"] == 10.452
    assert scan.scan_file(spirou)["mag_band"] == "H"
    write_tfits(nirps, instrument="NIRPS")          # no magnitude at all
    assert scan.scan_file(nirps)["mag"] is None
    assert scan.scan_file(nirps)["mag_band"] is None, \
        "a zero in a field nobody filled is not a magnitude"


def test_a_record_from_before_a_field_existed_is_read_again(tmp_path):
    """Adding a column must not throw a campaign's whole index away, and must
    not leave the new column empty for ever either."""
    root = root_with(tmp_path, PROXIMA=2)
    home = str(tmp_path / "home")
    index, _tally = scan.update(root, home=home)
    for record in index["objects"]["PROXIMA"]["files"].values():
        del record["mag"]                  # an index written before the column
    index, tally = scan.update(root, index=index, home=home)
    assert tally["read"] == 2 and tally["kept"] == 0
    assert all("mag" in f for f in index["objects"]["PROXIMA"]["files"].values())


def test_the_nights_two_campaigns_share_are_counted(tmp_path):
    """What decides whether a joint fit means anything: the atmosphere belongs
    to the night, so stars not observed together have no common sky to share."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    for name, days in (("GJ1", range(0, 20)), ("GJ3090", range(10, 30)),
                       ("PROXIMA", range(100, 120))):
        folder = tmp_path / name
        folder.mkdir()
        for i, day in enumerate(days):
            write_tfits(str(folder / ("%04dt.fits" % i)), mjd=60000.0 + day + 0.4)
    root = str(tmp_path)
    index, _tally = scan.update(root, home=str(tmp_path / "home"))

    assert len(scan.nights_of(index, "GJ1")) == 20
    common, counts = scan.shared_nights(index, ["GJ1", "GJ3090"])
    assert (len(common), counts) == (10, [20, 20]), "half of each"
    common, counts = scan.shared_nights(index, ["GJ1", "PROXIMA"])
    assert len(common) == 0, "different seasons, no common sky"
    common, _counts = scan.shared_nights(index, ["GJ1", "GJ3090", "PROXIMA"])
    assert len(common) == 0, "all three, so the one that shares nothing decides"
    assert scan.shared_nights(index, ["GJ1"]) == (set(), [20]), \
        "one object shares nothing with nobody"
    assert scan.shared_nights(index, ["NOBODY", "GJ1"]) == (set(), [20])


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


def test_the_barycentric_coverage_is_bins_that_hold_something(tmp_path):
    """The span and the coverage differ exactly where it matters: a campaign
    observed at two extremes and nowhere between has a wide span and almost no
    coverage, and it is the coverage that says whether the two frames can be
    told apart."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    for name, bervs in (("SPREAD", np.arange(-15.0, 15.1, 3.0)),
                        ("TWOENDS", np.array([-15.0, -14.0, 14.0, 15.0]))):
        folder = tmp_path / name
        folder.mkdir()
        for i, berv in enumerate(bervs):
            path = str(folder / ("%04dt.fits" % i))
            write_tfits(path, mjd=60000.0 + i)
            with fits.open(path, mode="update") as h:
                h[1].header["BERV"] = float(berv)
    root = str(tmp_path)
    index, _tally = scan.update(root, home=str(tmp_path / "home"))

    assert scan.bervs_of(index, "SPREAD").size == 11
    _edges, _counts, wide = scan.berv_coverage(index, ["SPREAD"])
    assert wide["span"] == 30.0
    assert wide["effective"] == 30.0, \
        "eleven points from -15 to 15 fill ten bins of 3 km/s"

    _edges, _counts, ends = scan.berv_coverage(index, ["TWOENDS"])
    assert ends["span"] == 30.0, "the same span"
    assert ends["effective"] == 6.0, "but two bins: the coverage says so"

    edges, counts, both = scan.berv_coverage(index, ["SPREAD", "TWOENDS"])
    assert set(counts) == {"SPREAD", "TWOENDS"}
    assert all(len(c) == len(edges) - 1 for c in counts.values()), "one grid"
    assert both["n"] == 15 and both["effective"] == 30.0
    assert both["objects"]["TWOENDS"]["effective"] == 6.0


def test_a_coverage_asked_of_nothing_is_zero_rather_than_a_crash(tmp_path):
    root = root_with(tmp_path, PROXIMA=2)     # written without a BERV keyword
    index, _tally = scan.update(root, home=str(tmp_path / "home"))
    edges, counts, summary = scan.berv_coverage(index, ["PROXIMA"])
    assert edges.size == 0 and counts == {} and summary["effective"] == 0.0
    assert scan.berv_coverage(index, [])[2]["span"] == 0.0


def test_the_ecliptic_latitude_says_what_a_target_can_ever_reach():
    """|BERV| <= 29.78 cos(beta): TOI-1452 at +80.5 deg can span 9.8 km/s and no
    more, whatever is observed. Checked against astropy on three targets."""
    for ra, dec, want in ((290.173917, 73.195, 80.49),
                          (217.376408, -62.67345, -44.77),
                          (133.781, 1.541, -15.21)):
        assert abs(scan._ecliptic_latitude(ra, dec) - want) < 0.02

    assert abs(scan.berv_limit(80.49) - 4.93) < 0.05, "TOI-1452, the hopeless one"
    assert abs(scan.berv_limit(0.0) - 29.78) < 0.01, "on the ecliptic, the most"
    assert abs(scan.berv_limit(90.0)) < 0.01, "at the pole, none at all"
    assert scan.berv_limit(None) is None


def test_the_coverage_says_what_the_sky_allowed_as_well(tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    folder = tmp_path / "POLE"
    folder.mkdir()
    for i, berv in enumerate((-2.0, 0.0, 2.0)):
        path = str(folder / ("%04dt.fits" % i))
        write_tfits(path, mjd=60000.0 + i)
        with fits.open(path, mode="update") as h:
            h[1].header["BERV"] = berv
            h[0].header["RA_DEG"] = 290.173917      # TOI-1452's, near the pole
            h[0].header["DEC_DEG"] = 73.195
    index, _tally = scan.update(str(tmp_path), home=str(tmp_path / "home"))
    _edges, _counts, summary = scan.berv_coverage(index, ["POLE"])
    assert summary["span"] == 4.0
    assert abs(summary["possible"] - 9.86) < 0.1, \
        "twice 29.78 cos(80.5 deg): what this target could ever have"
