"""
This is an attempt at a smarter implementation than directly following the setup from Knuth's 1974 article. The superpose step is replaced by an insert operation that inserts new closed sets into the family of current rank one at a time, superposing on the fly.
"""
function randomized_knuth_matroid_construction_v4(n, p, T=UInt16)::KnuthMatroid{T}
  r = 1
  pr = 0
  F = [Set(0)]
  E = 2^n - 1 # The set of all elements in E.

  while true
    to_insert = generate_covers_v2(F[r], n)

    # Apply coarsening to covers.
    if r <= length(p) && E ∉ to_insert # No need to coarsen if E is added.
      pr = p[r]
      while pr > 0
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

      for B in setdiff(F[r+1], A)
        if should_merge(A, B, F[r])
          push!(to_insert, A | B)
          setdiff!(F[r+1], [A, B])
          push!(F[r+1], A | B)
        end
      end
    end

    if E ∈ F[r+1]
      return KnuthMatroid{T}(n, F, [], Set(), Dict())
    end

    r += 1
  end
end