using DataFrames, CSV, DataFramesMeta

function makedists(munipop,municent)
    munij = leftjoin(munipop,municent, on = [:district,:municipality])
    tots = @by(munij, :district, :poptot = sum(:pop),:xtot = sum(:pop .* :centroid_x),:ytot = sum(:pop .* :centroid_y))
    data = @transform(tots, :distcode = :district, :pop = :poptot, :xcoord = :xtot ./ :poptot, :ycoord = :ytot ./ :poptot)
    @select!(data,:distcode,:pop,:xcoord,:ycoord)
end