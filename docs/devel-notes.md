## Notes

* https://google.github.io/styleguide/cppguide.html is taken _very roughly_ as
  a style guide
* `dnf install autoconf-archive` was required to get all autotools macros recognized.
* PFL "The Pitchfork" project layout is used
* .cxx is used for c++ source.
* .h is used for header files.
* external/xxHash/Makefile was edited to comment out the install and uninstall
  targets, and provide a no-op distdir target

## Release

* update version in src/version.h
* git commit src/version.h -m "..." && git push
* git tag [version]
* git push origin [version]
* ./autogen.sh  # AC_INIT picks up version from src/version.h
* make clean    # ensure no build artifacts in external/ which is included as-is in the tarball
* make dist     # creates streammd-[version].tar.gz
* In GitHub, draft new release from the tag and manually add the generated tarball
  to the release.

## Cross-compile for Apple ARM

* Clone [OSXcross](https://github.com/tpoechtrager/osxcross)
* Grab latest XCode from https://developer.apple.com/download/all/?q=xcode
  (you'll need an Apple ID)
* Follow the instructions [here](
  https://github.com/tpoechtrager/osxcross#packing-the-sdk-on-linux---method-1-xcode--80
) to package the SDK:
  
  and move the .tar.xz to `tarballs/` dir
* Set the SDK version (if > 1) and build the cross toolchain:
  ```
  SDK_VERSION=13.0 ./build.sh
  ```

* configure the build to use the cross toolchain
  ```
  cd path/to/streammd
  make clean
  CC=aarch64-apple-darwin22-clang CXX=aarch64-apple-darwin22-clang++ ./configure --host=aarch64-apple-darwin22
  ```

* xxHash external dependency is not an autotools project and autodetects the
  build platform from `uname` so we need to work around that:
  ```
  patch external/xxHash/Makefile external/xxHash/Makefile.osxcross.patch
  ```
* cross build the libxxhash.a first:
  ```
  cd external/xxHash
  make CC=arm64-apple-darwin22-clang CXX=arm64-apple-darwin22-clang++ AR=arm64-apple-darwin22-ar
  ```
* cross build the rest
  ```
  cd ../..
  make
  ```
* If all goes well you should end up with an arm64 executable
  ```
  file src/streammd
  # src/streammd: Mach-O 64-bit arm64 executable, flags:...
  ```
