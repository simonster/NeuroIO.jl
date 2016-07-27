using NeuralynxNCS, Base.Test

f = joinpath(dirname(@__FILE__), "TestFile.Ncs")
for ncs in (readncs(f), readncs(open(f)))
    @test ncs.header == "######## Neuralynx\r\nTest File\r\n"
    @test ncs.samples == [typemin(Int16):typemax(Int16);]
    dt = 1/32000
    @test ncs.times == [0:dt:dt*512*128-dt;]
    @test length(ncs.times.ranges) == 1
end
