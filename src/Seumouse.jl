module Seumouse

using DataFrames, CSVFiles, FileIO, SparseArrays

export mouse

function expressionframe(filename::AbstractString)::DataFrame
    f = open(filename)
    barcodevec = split(readline(f),",")
    barcodes = collect(b[2:end-1] for b in barcodevec)[2:end]
    cols = Pair{String,SparseVector{Float64,Int}}[]
    while !eof(f)
        l = split(readline(f),",")
        #first entry is our gene (surrounded by quotes)
        gene = l[1][2:end-1] |> String
        #the rest is the expression for each barcode
        sv = parse.(Float64,l[2:end]) |> sparse
        push!(cols,gene => sv)
    end
    #close the file
    close(f)
    #construct the dataframe
    DataFrame("barcode" => barcodes,cols...;copycols = false)
end

"""
```julia
mouse(meta,expression,other...)
```
Create a container for Seurat data from csv files. The required tables can be extracted
a Seurat object `obj` in R using expressions such as...
- meta: `as.data.frame(obj[[]])`
- expression: `as.data.frame(GetAssayData(object = obj[["integrated"]], layer = "data"))`
- umap: `as.data.frame(obj[["umap"]][[]])`
"""
function mouse(metafile::AbstractString,expressionfile::AbstractString,otherfiles::AbstractString...)
    nonexpressionframes = map([metafile,otherfiles...]) do f
        frame = load(f) |> DataFrame
        rename!(frame,1 => :barcode)
    end
    df = innerjoin(nonexpressionframes...,on = :barcode)
    
    #get our expression object
    ef = expressionframe(expressionfile)

    #now we want to make sure that the barcodes of df are in the same order as ef, because we want to
    #hcat without copying any of the rows
    #this little helper dict will give the position of the barcode in ef.barcode so we can sort the other frame
    barcoderank = Dict(b => i for (i,b) in enumerate(ef.barcode))

    #this line gives the permutation which will put df in the same barcode order as ef
    p = sortperm(df.barcode,by = (b) -> barcoderank[b])

    #reorder df if it isn't already ordered
    df = df[p,:]

    @assert all(df.barcode .== ef.barcode) "barcode columns do not match"

    #hcat without copying and return. need to strip the barcode column from ef so all cols are unique
    hcat(df,select(ef,Not(:barcode)),copycols=false)
end

end # module Seumouse
