# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Configure (out-of-source build required)
cmake -B builds -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build builds --config Release

# Run all tests
cd builds && ctest -C Release

# Run GoogleTest binary directly (for verbose output)
./builds/test/Tiled_Matrix_gtest

# Generate Doxygen docs (if Doxygen is installed)
cmake --build builds --target doc_doxygen
```

> In-source builds are blocked by CMakeLists.txt. Always build in a separate directory (e.g., `builds/` or `cmake-build-debug/`).

## Architecture

This is **Tiled_Matrix**, a C++20 header/library for tiled (blocked) matrix storage.

### Core design: tiled storage layout

A `Matrix<T>` stores an `m × n` matrix divided into `p × q` tiles of size `mb × nb`. All data is held in a single flat `std::unique_ptr<T[]>` (`top_`), with element access converting 2D or tile-local coordinates into a 1D index via `convertTileToArray()`.

The `Ordering` enum controls two independent layout choices:
- **Tile matrix order**: row-major or column-major arrangement of tiles within the full matrix
- **Tile internal order**: row-major or column-major layout of elements within each tile

Default ordering is `TileMatrixColumnMajorTileColumnMajor`.

### Key files

| File | Role |
|------|------|
| `src/matrix.hpp` | Full class definition, constructors, inline getters, `operator<<` |
| `src/matrix.cpp` | Method implementations; explicitly instantiated for `float` and `double` only |

### Element access

- `A(i, j)` — logical row/col access (converts to tile coords internally)
- `A[k]` — flat 1D index into the underlying array
- `A.elm(ti, tj)` — pointer to start of tile `(ti, tj)`
- `A.elm(ti, tj, i, j)` — pointer to element `(i,j)` within tile `(ti,tj)`

### Boundary tiles

`mb(ti, tj)` and `nb(ti, tj)` (the two-argument overloads) return the actual height/width of a tile, accounting for remainder rows/columns in the last tile. The single-argument getters `mb()` / `nb()` return the nominal tile size.

### Template instantiation

`matrix.cpp` ends with explicit instantiations for `float` and `double`. To support other numeric types (e.g., `std::complex<double>`), add an explicit instantiation there.

### Tests

`test/matrix_googletest.cpp` — GoogleTest suite covering `equal`, `plus`, `minus`, `elm`, and `zero` for both `float` and `double`.  
`test/Matrix_test.cpp` — standalone executable (no test framework).

GoogleTest is fetched via CMake `FetchContent` at configure time (tag `v1.13.0`).
