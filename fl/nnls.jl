using NPZ
import NonNegLeastSquares

shot = parse(Int, ARGS[1])
segment_count = parse(Int, ARGS[2])

for segment=0:segment_count-1
    # Load data and count dimensions
    println("STATUS: Began working on segment $segment")
    fl_matrix = npzread("../cache/fl_matrix_Xpt_$(shot)_$(segment).npy")
    frames = npzread("../cache/frames_Xpt_$(shot)_$(segment).npy")
    fl_count = size(fl_matrix)[2]
    frame_count = size(frames)[1]

    # Distribute NNLS reconstructions among workers
    fl_emissivities = pmap((a1, a2)->NonNegLeastSquares.nonneg_lsq(a1, a2; alg=:admm), [fl_matrix for i=1:frame_count], [frames[j, :] for j=1:frame_count])

    # Convert abstract-typed array to matrix for npzwrite
    fl_emissivities_array = Array{Float64}(frame_count, fl_count)
    for i=1:frame_count
        for j=1:fl_count
            fl_emissivities_array[i, j] = fl_emissivities[i][j]
        end
    end

    padded_segment = dec(segment, 2)
    npzwrite("../cache/fl_emissivities_Xpt_$(shot)_$(padded_segment).npy", Array(fl_emissivities_array))
end
