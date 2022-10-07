## Notes

* `dnf install autoconf-archive` was required to get all autotools macros recognized.
* PFL "The Pitchfork" project layout is used
* .cxx is used for c++ source.
* .h is used for header files.

## Release

* update version in src/version.h
* git commit src/version.h -m "..." && git push
* git tag [version]
* git push origin [version]
* ./autogen.sh  # AC_INIT picks up version from src/version.h
* make dist     # creates streammd-[version].tar.gz
* In GitHub, draft new release from the tag and manually add the generated tarball
  to the release.
