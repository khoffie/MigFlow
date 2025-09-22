getageyear(r) = r.mdl.data.age, r.mdl.data.year

function grid_position(k::Int, ncols::Int=2)
    i = div(k - 1, ncols) + 1   # row index
    j = mod(k - 1, ncols) + 1   # column index
    return (i, j)
end

function prettytitle!(f, title)
    titlelayout = GridLayout(f[0, 1], tellwidth = false)
    Label(titlelayout[1, 1], title, fontsize = 15, font = "TeX Gyre Heros Bold Makie")
    rowgap!(titlelayout, 0)
    rowgap!(f.layout, 0)
    colgap!(f.layout, -50)
end

function hideall!(ax)
    tightlimits!(ax)
    hidedecorations!(ax)
    hidespines!(ax)
end
