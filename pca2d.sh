#!/bin/sh
# One command to be in the right place: leave whatever conda environment is
# active, enter pca2d-preclean, and say what there is to run.
#
# It has to be SOURCED, since a child process cannot change its parent's
# environment. Put this line in ~/.zshrc (or ~/.bashrc), with the path of this
# checkout:
#
#     alias pca2d='source /path/to/pca2d_preclean/pca2d.sh'
#
# Then `pca2d` from anywhere.
#
# The environment holds BOTH codes, this one and LBL, which is why leaving the
# one that happens to be active matters: having to deactivate one to run the
# other is how a t.fits gets measured by the wrong version of something.

ENV_NAME=${PCA2D_ENV:-pca2d-preclean}

# --- this file must be sourced -----------------------------------------------
_pca2d_sourced=0
if [ -n "${ZSH_VERSION:-}" ]; then
    case ${ZSH_EVAL_CONTEXT:-} in *:file*) _pca2d_sourced=1 ;; esac
    _pca2d_self=${(%):-%N}
elif [ -n "${BASH_VERSION:-}" ]; then
    [ "${BASH_SOURCE[0]}" != "${0}" ] && _pca2d_sourced=1
    _pca2d_self=${BASH_SOURCE[0]}
else
    _pca2d_sourced=1            # unknown shell: assume the alias is right
    _pca2d_self=$0
fi

# --- where this checkout is --------------------------------------------------
# resolved from the file being sourced, so the alias may name any checkout, and
# absolute, so the line it suggests for an rc file is one that works from
# anywhere
case ${_pca2d_self} in
    /*) PCA2D_HOME=$(cd "$(dirname "${_pca2d_self}")" && pwd) ;;
    *)  PCA2D_HOME=$(cd "$(dirname "${PWD}/${_pca2d_self}")" && pwd) ;;
esac
_pca2d_path="${PCA2D_HOME}/$(basename "${_pca2d_self}")"

if [ "$_pca2d_sourced" -eq 0 ]; then
    printf '%s\n' "source this file rather than running it, or a child process"
    printf '%s\n' "changes its own environment and exits. The alias to keep:"
    printf '\n    alias pca2d='"'"'source %s'"'"'\n\n' "${_pca2d_path}"
    printf '%s\n' "and to write it once:"
    printf "    echo \"alias pca2d='source %s'\" >> ~/.zshrc\n" "${_pca2d_path}"
    exit 1
fi

# --- colours, and only on a terminal ----------------------------------------
if [ -t 1 ]; then
    _c_head=$(printf '\033[1;34m'); _c_cmd=$(printf '\033[1;32m')
    _c_dim=$(printf '\033[2m');     _c_warn=$(printf '\033[33m')
    _c_err=$(printf '\033[31m');    _c_off=$(printf '\033[0m')
else
    _c_head=; _c_cmd=; _c_dim=; _c_warn=; _c_err=; _c_off=
fi

# --- conda, however deep the shell already is --------------------------------
# `conda activate` needs conda's SHELL FUNCTION, not the executable on PATH: a
# shell that has only the executable answers "Run 'conda init' first". A path
# back from command -v means that is all there is, so load the function.
_conda=$(command -v conda 2>/dev/null)
case ${_conda} in
    ""|/*)
        _base=${CONDA_EXE:-${_conda}}
        [ -n "${_base}" ] && \
            . "$(dirname "$(dirname "${_base}")")/etc/profile.d/conda.sh" 2>/dev/null
        unset _base
        ;;
esac
unset _conda

if ! command -v conda >/dev/null 2>&1; then
    printf '%sno conda in this shell%s: install it, or run `conda init` once.\n' \
        "${_c_err}" "${_c_off}"
else
    # leave every environment that is active, nested ones included, then enter
    # this one. Deactivating rather than activating on top of whatever was
    # there keeps a second environment's paths out of this one.
    _left=
    while [ "${CONDA_SHLVL:-0}" -gt 0 ]; do
        _here=${CONDA_DEFAULT_ENV:-base}
        # named once, however many times a shell stacked it
        case " ${_left} " in
            *" ${_here} "*) : ;;
            *) [ "${_here}" != "${ENV_NAME}" ] && _left="${_left} ${_here}" ;;
        esac
        conda deactivate || break
    done
    unset _here
    if conda activate "${ENV_NAME}" 2>/dev/null; then
        [ -n "${_left}" ] && printf '%sleft%s%s\n' \
            "${_c_dim}" "${_left}" "${_c_off}"
    else
        printf '%sno environment called %s%s. Make it once, from this folder:\n' \
            "${_c_err}" "${ENV_NAME}" "${_c_off}"
        printf '    %sconda env create -f %s/environment.yml%s\n' \
            "${_c_cmd}" "${PCA2D_HOME}" "${_c_off}"
    fi
    unset _left
fi

# --- what there is to run ----------------------------------------------------
printf '\n%spca2d-preclean%s  %stwo frames, one fit, then LBL%s\n' \
    "${_c_head}" "${_c_off}" "${_c_dim}" "${_c_off}"
printf '%s  environment %s%s   %spython %s%s\n' \
    "${_c_dim}" "${_c_off}" "${CONDA_DEFAULT_ENV:-none}" \
    "${_c_dim}" "$(python -c 'import platform;print(platform.python_version())' 2>/dev/null || echo '?')" "${_c_off}"
printf '%s  checkout    %s%s\n' "${_c_dim}" "${_c_off}" "${PCA2D_HOME}"
printf '%s  config      %s%s/config.yaml\n\n' "${_c_dim}" "${_c_off}" "${PCA2D_HOME}"

printf '  %spca2d-gui%s                                   the window: pick targets, run, watch it\n' \
    "${_c_cmd}" "${_c_off}"
printf '  %spca2d-preclean --object TOI-2120%s            one target, t.fits to velocities\n' \
    "${_c_cmd}" "${_c_off}"
printf '  %spca2d-preclean --objects A,B,C --n-star 0%s   several stars, one observer basis\n' \
    "${_c_cmd}" "${_c_off}"
printf '  %spca2d-preclean --object X --dry-run%s         resolve everything, touch nothing\n' \
    "${_c_cmd}" "${_c_off}"
printf '  %s./check.sh%s                                  the regression net, 13 s, in the checkout\n\n' \
    "${_c_cmd}" "${_c_off}"

printf '%s  spectra go in   <input root>/<object>/, which is never written to\n' "${_c_dim}"
printf '  everything else lands under the output root, one folder per run\n'
printf '  to read first   docs/options.md, the nominal and what each choice was worth%s\n\n' \
    "${_c_off}"

unset _pca2d_sourced _pca2d_self _pca2d_path _c_head _c_cmd _c_dim _c_warn _c_err _c_off
