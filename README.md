# map2ava
Convert read-to-ref mapping to read-vs-read overlapping. If the mapping regions of two reads overlap, then the overlapping region/positions of these two reads can be calculated by their mapping CIGARs.

## Installation

#### Dependencies

`map2ava` is written in rust, try below commands (no root required) or see [here](https://www.rust-lang.org/tools/install) to install `Rust` first.
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

#### Download and install

```
git clone https://github.com/moold/map2ava.git
cd map2ava && cargo build --release
```

## Usage
`./target/release/map2ava --thread 20 input.map.bam > input.overlap.paf`

## Parameters
Use `./target/release/map2ava -h` to see options.