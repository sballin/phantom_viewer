using NPZ
import NonNegLeastSquares

shot = parse(Int, ARGS[1])
segment_count = parse(Int, ARGS[2])
smoothing_param = parse(Int, ARGS[3])

for segment=0:segment_count-1
    # Load data and count dimensions
    println("STATUS: Began working on segment $segment")
    fl_matrix = npzread("../cache/fl_matrix_Xpt_$(shot)_$(segment).npy")
    frames = npzread("../cache/frames_Xpt_$(shot)_$(segment).npy")
    
    fl_count = size(fl_matrix)[2]
    frame_count = size(frames)[1]
    if smoothing_param == 0
        method = :fnnls
    else
        method = :admm
    end

    # Distribute NNLS reconstruction jobs among workers
    tic()
    fl_emissivities = pmap((a1, a2)->NonNegLeastSquares.nonneg_lsq(a1, a2; alg=method), [fl_matrix for i=1:frame_count], [frames[j, :] for j=1:frame_count])
    toc()

    # Convert abstract-typed array to matrix for npzwrite
    fl_emissivities_array = Array{Float64}(frame_count, fl_count)
    for i=1:frame_count
        for j=1:fl_count
            fl_emissivities_array[i, j] = fl_emissivities[i][j]
        end
    end
    println(sum(fl_emissivities_array))

    padded_segment = dec(segment, 2)
    npzwrite("../cache/fl_emissivities_Xpt_$(shot)_sp$(smoothing_param)_$(padded_segment).npy", Array(fl_emissivities_array))
end
