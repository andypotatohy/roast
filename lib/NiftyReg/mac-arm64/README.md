# NiftyReg macOS Apple Silicon runtime

These binaries come from the official KCL-BMEIS NiftyReg v2.0.0 release asset:

https://github.com/KCL-BMEIS/niftyreg/releases/download/v2.0.0/NiftyReg-macOS-v2.0.0.zip

Release checksum:

```text
sha256 83d75c72a0d3b61102beaddc779d7b69323a3fdcfc30d65b7e8c0b08ffb37702
```

The binaries are arm64 Mach-O executables for Apple Silicon.

Runtime dylibs needed by this release are bundled in `lib/` and the executable load paths have been rewritten to use `@executable_path/lib`, so users do not need Homebrew to run these NiftyReg binaries.

Bundled runtime libraries:

- `libgomp.1.dylib`
- `libstdc++.6.dylib`
- `libgcc_s.1.1.dylib`
- `libpng16.16.dylib`

NiftyReg is distributed under the BSD-3-Clause license; see `LICENSE.txt` in this directory. Runtime library licenses are included under `licenses/`.
