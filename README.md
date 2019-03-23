A simple Julia type for holding Fisher matrices along with parameter names so
that they can be added together and the rows and columns are automatically
aligned. 

To create a Fisher matrix:

```julia
julia> f = FisherMatrix([1 0; 0 2],["param1","param2"])
Fisher.FisherMatrix{Int64,ASCIIString}, std-dev of 2 parameters: 
  param2 => 0.7071067811865476 
  param1 => 1.0
```

The array itself and the names are stored under `f.fish` and `f.names`. We can
create a second  fisher matrix and add them together. Note the rows/columns are
automatically aligned. 

```julia
julia> f2 = FisherMatrix([3 0; 0 4],["param2","param3"])
Fisher.FisherMatrix{Int64,ASCIIString}, std-dev of 2 parameters: 
  param2 => 0.5773502691896257 
  param3 => 0.5 
  
julia> f + f2
Fisher.FisherMatrix{Float64,ASCIIString}, std-dev of 3 parameters: 
  param2 => 0.4472135954999579 
  param1 => 1.0 
  param3 => 0.5 
```

The standard deviations (the sqrt's of the diagonal of the inverse) can be
extracted with the `stds` method
```julia
julia> stds(f + f2)
Dict{ASCIIString,Float64} with 3 entries:
  "param2" => 0.4472135954999579
  "param1" => 1.0
  "param3" => 0.5
```
