using NPZ
import NonNegLeastSquares

shot = parse(Int, ARGS[1])
segment_count = parse(Int, ARGS[2])

for segment=0:segment_count-1
    println("STATUS: Began working on segment $segment")
    fl_matrix = npzread("../cache/fl_matrix_Xpt_$(shot)_$(segment).npy")
    frames = npzread("../cache/frames_Xpt_$(shot)_$(segment).npy")
    fl_count = size(fl_matrix)[2]
    frame_count = size(frames)[1]

    mt = fl_matrix'
    mtm = mt * fl_matrix
    fl_emissivities = SharedArray(Float64, (frame_count, fl_count))
    @sync @parallel for frame_i=1:frame_count
        fl_emissivities[frame_i, :] = NonNegLeastSquares.fnnls(mtm, mt * frames[frame_i, :])
            println("STATUS: NNLS is at frame $frame_i out of $frame_count")
    end

    npzwrite("../cache/fl_emissivities_Xpt_$(shot)_$(segment).npy", Array(fl_emissivities))
end
