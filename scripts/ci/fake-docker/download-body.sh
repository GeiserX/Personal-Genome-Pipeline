#!/usr/bin/env bash
# download-body.sh: what the fake curl and wget "download". Sourced by both.
#
# fake_body URL prints the body for URL:
#   - the file ${FAKE_DOWNLOAD_DIR}/<base name of URL>, when a case put one
#     there (a real tarball, a checksum file with a wrong value...);
#   - for a URL ending in .md5, the md5 of the body of the URL without .md5,
#     so a script that verifies a download against the published md5 passes;
#   - a gzip-compressed placeholder for .gz, .bgz and .tgz, so `gzip -t` passes;
#   - a one-line placeholder for everything else.

fake_md5() {
  if command -v md5sum >/dev/null 2>&1; then md5sum | awk '{print $1}'
  elif command -v md5 >/dev/null 2>&1; then md5 -q
  else openssl md5 -r | awk '{print $1}'; fi
}

fake_body() {
  local url=$1 name
  name=$(basename "$url")
  if [ -n "${FAKE_DOWNLOAD_DIR:-}" ] && [ -f "${FAKE_DOWNLOAD_DIR}/${name}" ]; then
    cat "${FAKE_DOWNLOAD_DIR}/${name}"
    return
  fi
  case "$url" in
    *.md5) printf '%s  %s\n' "$(fake_body "${url%.md5}" | fake_md5)" "${name%.md5}" ;;
    *.gz|*.bgz|*.tgz) printf 'placeholder for %s\n' "$url" | gzip -nc ;;
    *) printf 'placeholder for %s\n' "$url" ;;
  esac
}
