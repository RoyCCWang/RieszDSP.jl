import PythonPlot as PLT
using Serialization

PLT.close("all")
fig_num = 1

m = 1
ts_tar, ys_tar, ts_ref, ys_ref = deserialize(joinpath("tmp", "xic_$(m)"))

PLT.figure(fig_num)
fig_num += 1
PLT.plot(ts_ref, ys_ref, "x", label="reference")
PLT.plot(ts_tar, ys_tar, "x", label="target")
PLT.legend()
PLT.title("XIC $m data")

# The levels were itnroduced by itp extrapolation.

nothing
