"""
Sixth implementation of random-KMC, in which a rank table is used to keep track of set ranks, and the covers and enlargements are added one at a time, ensuring the matroid properties at all times.
"""
function random_kmc_v6(n, p, T=UInt16)::KnuthMatroid{T}
  r = 1
  pr = 0
  F::Vector{Set{T}} = [Set(T(0))]
  E::T = BigInt(2)^n-1
  rank = Dict{T, UInt8}(0=>0)

  while E ∉ F[r]
    # Create empty set.
    push!(F, Set())
    
    # Generate minimal closed sets for rank r+1.
    for y in F[r] # y is a closed set of rank r.
      t = E - y # The set of elements not in y.
      # Find all sets in F[r+1] that already contain y and remove excess elements from t.
      for x in F[r+1]
        if (x & y == y) t &= ~x end
      end
      # Insert y ∪ a for all a ∈ t.
      while t > 0
        x = y|(t&-t)
        add_set!(x, F, r, rank)
        t &= ~x
      end
    end

    if E ∈ F[r+1]
      break
    end

    if r <= length(p)
      # Apply coarsening.
      pr = p[r]
      while pr > 0 && E ∉ F[r+1]
        A = rand(F[r+1])
        t = E-A
        one_element_added::Vector{T} = []
        while t > 0
          x = A|(t&-t)
          push!(one_element_added, x)
          t &= ~x
        end
        Acupa = rand(one_element_added)
        setdiff!(F[r+1], A)
        add_set!(Acupa, F, r, rank)
        pr -= 1
      end
    end

    r += 1
  end

  return KnuthMatroid{T}(n, F, [], Set(), rank)
end

function add_set!(x, F, r, rank)
  if x in F[r+1] return end
  for y in F[r+1]
    if haskey(rank, x&y) && rank[x&y]<r
      continue
    end

    # x ∩ y has rank > r, replace with x ∪ y.
    setdiff!(F[r+1], y)
    return add_set!(x|y, F, r, rank)
  end

  push!(F[r+1], x)
  rank[x] = r
end