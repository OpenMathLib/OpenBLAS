# Meson build

OpenBLAS also offers Meson build tool support. Currently `Linux` and `Macos` systems are supported.
In terms of supported architectures please refer to the table below:

|        | subarchitectures                    |
|--------|-------------------------------------|
| x86_64 | haswell, sandybridge, skylakex, zen |
| armv8  | armv8                               |

To use Meson as a build tool, clone the repo, make sure that `meson` command is available (it can
be installed within a new conda environment e.g. `openblas-dev`):

```bash
meson setup build --buildtype release

meson compile -C build

meson test -C build -v
```

In case any of the `meson.build` were changed, use `--reconfigure` to regenerate targets:

```bash
meson setup build --reconfigure --buildtype release
```

## Implementation details

Meson build aims to replicate Makefile setup as much as possible, while keeping the `meson.build`
files well-structured and straightforward. In general each compiled object is associated with
a **dict** entry (in `meson.build` files inside `kernel`, `driver`, and `interface` directories).

Inside `kernel` directory each entry is specified by a function name, data type, and (possibly)
extension and points to the source file that should be used for building the object.

Both `kernel`'s Makefiles and `meson.build` files have hierarchical structure, but in Meson build
they're purely code-driven. Assuming we have a `HASWELL` machine the search order for the source
file will be:

```
1. x86_64_haswell_dict

2. x86_64_base_dict

3. base_dict
```

This allows to easly override generic configuration with architecture-specific details.
In Makefiles the overriding logic implementation might appear cluttered, with multiple `if`
statements and architecture-specific override rules.

With Meson the search order is established by `search_order` variable which allows to quickly
determine the input file for any function.
Object-specific build flags are stored separately, defined at an extension level, inside
`meson.build` and `kernel/meson.build`.
