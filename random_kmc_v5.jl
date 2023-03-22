"""
Fifth implementation of random-KMC. This one uses a dictionary to keep track of previously seen sets.
"""
function random_kmc_v5(n, p, T=UInt16)::KnuthMatroid{T}
  r = 1
  pr = 0
  F::Vector{Set{T}} = [Set(T(0))]
  E = 2^n-1
  rank = Dict{T, UInt8}(0=>0) # The rank table maps from the representation of a set to its assigned rank.

  while true
    to_insert = generate_covers_v2(F[r], n)

    # Apply coarsening to covers.
    if r <= length(p)
      pr = p[r]
      while length(to_insert) > 0 && pr > 0 && E ∉ to_insert # No need to coarsen if E is added.
        A = rand(to_insert)
        a = random_element(E - A)
        to_insert = setdiff(to_insert, A) ∪ [A | a]
        pr -= 1
      end
    end

    # Superpose.
    push!(F, Set()) # Add F[r+1].
    while length(to_insert) > 0
      A = pop!(to_insert)
      push!(F[r+1], A)
      rank[A] = r

      for B in setdiff(F[r+1], A)
        if !haskey(rank, A&B) || rank[A&B] >= r
          # Update insert queue.
          push!(to_insert, A | B)

          # Update F[r+1].
          setdiff!(F[r+1], [A, B])
          push!(F[r+1], A | B)

          # Update rank table.
          rank[A|B] = r
          break
        end
      end
    end

    if E ∈ F[r+1]
      return KnuthMatroid{T}(n, F, [], Set(), rank)
    end

    r += 1
  end
end