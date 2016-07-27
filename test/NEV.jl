using NeuroIO, MAT, BinDeps, Base.Test

# Unfortunately, these test files come from Black Rock Microsystems and
# are assumed not to be redistributable
testdir = joinpath(dirname(@__FILE__), "NSS Data")
if !isdir(testdir)
    parentdir = dirname(@__FILE__)
    bundlezip = joinpath(parentdir, "NSS Data.zip")
    run(download_cmd("http://cl.ly/1B3y1C0U380G/download/NSS%20Data.zip", bundlezip))
    cd(parentdir) do
        run(`unzip $bundlezip`)
    end
end

if !isfile(joinpath(testdir, "nev.mat"))
    using MATLAB
    npmkdir = joinpath(dirname(@__FILE__), "NPMK")
    if !isdir(npmkdir)
        cd(dirname(@__FILE__)) do
            run(`git clone https://github.com/BlackrockMicrosystems/NPMK.git`)
        end
    end
    @eval mat"""
    addpath(genpath($(joinpath(npmkdir, "NPMK"))))
    cd($testdir);
    NEV = openNEV('./SampleNSSData.nev');
    save('nev.mat', 'NEV');
    NS5 = openNSx('./SampleNSSData.ns5');
    save('ns5.mat', 'NS5');
    """
end

nev = NeuroIO.NEV.readnev(open("NSS Data/SampleNSSData.nev"))
nev_mat = matread("NSS Data/nev.mat")["NEV"]
@test nev.header.bytes_in_headers == nev_mat["MetaTags"]["HeaderOffset"]
@test nev.header.bytes_in_data_packets == nev_mat["MetaTags"]["PacketBytes"]
@test nev.header.time_resolution_of_time_stamps == nev_mat["MetaTags"]["TimeRes"]
@test nev.header.time_resolution_of_samples == nev_mat["MetaTags"]["SampleRes"]
@test nev.header.application_to_create_file == nev_mat["MetaTags"]["Application"]

info = nev_mat["ElectrodesInfo"]
spikes = nev_mat["Data"]["Spikes"]
for i = 1:length(info["ElectrodeID"])
    id = info["ElectrodeID"][i]
    ch = nev.spike_channels[id]
    @test ch.physical_connector == info["ConnectorBank"][i][1]
    @test ch.connector_pin == info["ConnectorPin"][i]
    # What is this?
    # @test ch.digitization_factor == info["DigitalFactor"][i]
    @test ch.energy_threshold == info["EnergyThreshold"][i]
    @test ch.high_threshold == info["HighThreshold"][i]
    @test ch.low_threshold == info["LowThreshold"][i]
    @test ch.number_of_sorted_units == info["Units"][i]
    @test ch.bytes_per_waveform == info["WaveformBytes"][i]
    @test ch.label == rstrip(string(info["ElectrodeLabel"][i]...), '\0')
    @test ch.high_pass.cutoff == info["HighFreqCorner"][i]
    @test ch.high_pass.order == info["HighFreqOrder"][i]
    @test ch.high_pass.filter_type == info["HighFilterType"][i]
    @test ch.low_pass.cutoff == info["LowFreqCorner"][i]
    @test ch.low_pass.order == info["LowFreqOrder"][i]
    @test ch.low_pass.filter_type == info["LowFilterType"][i]

    spike_index = find(spikes["Electrode"] .== id)
    @test ch.times == spikes["TimeStamp"][spike_index]/nev.header.time_resolution_of_time_stamps
    @test ch.unit_numbers == spikes["Unit"][spike_index]
    @test ch.waveforms == spikes["Waveform"][:, spike_index]
end

ns5 = NeuroIO.NEV.readnsx(open("NSS Data/SampleNSSData.ns5"))
ns5_mat = matread("NSS Data/ns5.mat")["NS5"]
@test ns5.header.label == rstrip(string(ns5_mat["MetaTags"]["SamplingLabel"]...), '\0')
@test ns5.header.period == 30000/ns5_mat["MetaTags"]["SamplingFreq"]
@test ns5.header.time_resolution_of_time_stamps == ns5_mat["MetaTags"]["TimeRes"]

info = ns5_mat["ElectrodesInfo"]
for i = 1:length(info["ElectrodeID"])
    id = info["ElectrodeID"][i]
    ch = ns5.continuous_channels[id]
    @test ch.label == rstrip(string(info["Label"][i]...), '\0')
    @test ch.physical_connector == info["ConnectorBank"][i][1]
    @test ch.connector_pin == info["ConnectorPin"][i]
    @test ch.min_digital_value == info["MinDigiValue"][i]
    @test ch.max_digital_value == info["MaxDigiValue"][i]
    @test ch.min_analog_value == info["MinAnalogValue"][i]
    @test ch.max_analog_value == info["MaxAnalogValue"][i]
    @test ch.units == rstrip(string(info["AnalogUnits"][i]...), '\0')
    @test ch.high_pass.cutoff == info["HighFreqCorner"][i]
    @test ch.high_pass.order == info["HighFreqOrder"][i]
    @test ch.high_pass.filter_type == info["HighFilterType"][i]
    @test ch.low_pass.cutoff == info["LowFreqCorner"][i]
    @test ch.low_pass.order == info["LowFreqOrder"][i]
    @test ch.low_pass.filter_type == info["LowFilterType"][i]

    @test ch.samples == vec(ns5_mat["Data"][i, :])
end
