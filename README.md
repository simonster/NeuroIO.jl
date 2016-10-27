# NeuroIO

[![Build Status](https://travis-ci.org/simonster/NeuroIO.jl.svg?branch=master)](https://travis-ci.org/simonster/NeuroIO.jl)

## Quick Start

```julia
using NeuroIO
f = open("/path/to/file")
x = readplx(f)                          # or readnev, readnsx

NeuroIO.times(x.spike_channels[1])      # spike times for channel 1
unitnumbers(x.spike_channels[1])        # unit numbers for spikes
NeuroIO.data(x.spike_channels[1])       # waveforms
NeuroIO.label(x.spike_channels[1])      # channel label
voltagemultiplier(x.spike_channels[1])  # voltage multiplier

NeuroIO.data(x.continuous_channels[1])  # continuous data for channel 1
NeuroIO.times(x.continuous_channels[1]) # sample times for channel 1

NeuroIO.data(x.event_channels[257])     # event codes for channel 257
NeuroIO.times(x.event_channels[257])    # event times

validchannels(x.continuous_channels)    # recorded continuous channels
```

## TODO

- Write NEV/NSX and NCS files
- Improve test coverage
- Improve `show`
- Documentation
