"""
Repack  GeGeGe zipped, chunks 1e7 to unzipped, chunks 10000
This speeds up operation time significantly (factor of 3) while
file size is about 2-times bigger

"""

using HDF5
using Printf
using Dates

function process(din::HDF5Dataset, dout::HDF5Dataset, 
                 n::Int64, o_chunk_size::Int64)
    left_pos = 1
    is_something_left = true
    is_something_left = true
    i_chunk_size = HDF5.get_chunk(din)[2]
    n_chunks = div(i_chunk_size, o_chunk_size)

    t0 = Dates.Time(Dates.now())

    while is_something_left
        right_pos = left_pos + i_chunk_size - 1
        if right_pos > n
            right_pos = n
            is_something_left = false
        elseif left_pos >= right_pos
            break
        end
        dout[:, left_pos:right_pos] = din[:, left_pos:right_pos]
        left_pos = right_pos + 1

        t1 = Dates.Time(Dates.now())
        dt = t1 - t0
        print("\r", round(left_pos / n * 100, digits=1), "% ",
              round(dt.value * 1e-9, digits=1), " s (",
              round(dt.value * 1e-9 * n / left_pos, digits=0), " s)")
    end
    println()
end

fin = h5open("u238_m_ggg.h5", "r")
fout = h5open("new.h5", "w")

gin = fin["GeGeGe"]
gout = create_group(fout, "GeGeGe")

chunk_size = 10_000

for i in 3:9
    s = @Printf.sprintf("%d", i)
    println("m = ", s)
    din = gin[s]
    n = size(din)[2]
    dout = create_dataset(gout, s, datatype(UInt32), dataspace(7, n),
                       chunk=(7, chunk_size))
    process(din, dout, n, chunk_size)
end
close(fin)
close(fout)
