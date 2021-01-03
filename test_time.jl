"""
Test timing distrubutions 
* gamma-time (for all events broken into multiplicities
* gamma-time for 606-423 gate (136Te)

"""
import HDF5

include("Nucubes.jl")
import .Nucubes

fin_u = HDF5.h5open("u238_m_ggg.h5", "r")
fin_th = HDF5.h5open("th232_m_ggg.h5", "r")
fout = HDF5.h5open("time_test.h5", "w")

for m in 3:9
    gate_m = [m, m]

    tg = Nucubes.TGate([2048, 512], 100, 1000)
    matrix = Nucubes.process(fin_u, gate_m, tg, Nucubes.select_t!)
    println("u ", m, " total counts ", sum(matrix))
    d = HDF5.create_dataset(fout, "u_dt_E_$m", HDF5.datatype(Int64), 
                            HDF5.dataspace(size(matrix)))
    d[:, :] = matrix[:, :]

    ggt = Nucubes.GGTGate([2048, 512], [605.6, 607.6], [422.4, 424.4], 
                          100, 1000)
    matrix = Nucubes.process(fin_u, gate_m, ggt, Nucubes.select_ggt!)
    println("u 606-423 ", m, " total counts ", sum(matrix))
    d = HDF5.create_dataset(fout, "u_606_423_E_t_m$m", HDF5.datatype(Int64), 
                            HDF5.dataspace(size(matrix)))
    d[:, :] = matrix[:, :]

    matrix = Nucubes.process(fin_th, gate_m, tg, Nucubes.select_t!)
    println("th ", m, " total counts ", sum(matrix))
    d = HDF5.create_dataset(fout, "th_dt_E_$m", HDF5.datatype(Int64), 
                            HDF5.dataspace(size(matrix)))
    d[:, :] = matrix[:, :]

    matrix = Nucubes.process(fin_th, gate_m, ggt, Nucubes.select_ggt!)
    println("606-423 th ", m, " total counts ", sum(matrix))
    d = HDF5.create_dataset(fout, "th_606_423_E_t_m$m", HDF5.datatype(Int64), 
                            HDF5.dataspace(size(matrix)))
    d[:, :] = matrix[:, :]
end
close(fin_u)
close(fin_th)
close(fout)
