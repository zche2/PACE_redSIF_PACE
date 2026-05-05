using JLD2
using Plots

# load tropomi retrieval results
svd_file = "/kiwi-data/Data/groupMembers/zhe2/MyProjects/TROPOMI_SIF/svd/svd_5_20.0_200.0_2.0_80.0_b6__v10.jld2"
jldopen(svd_file, "r") do f
    println(keys(f))
end
# load svd single stored object
@load svd_file single_stored_object
fieldnames(typeof(single_stored_object))

single_stored_object.dates
single_stored_object.svd_file_info
single_stored_object.wavelength_index_bounds
cross_track_num, radiance_bin, 
size(single_stored_object.svd)

# access the first svd cell
U, S, Vt = single_stored_object.svd[1, 1, 1]

println("Size of U: ", size(U))
println("Length of S: ", length(S))
println("Size of Vt: ", size(Vt))