"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the size of the universe. Bit-based.
"""
function generate_covers_v2(F_r, n)
  Set([A | 1 << i for A ∈ F_r for i in 0:n-1 if A & 1 << i === 0])
end

"""
Returns whether the intersection of A and B is contained within
some C in F_prev.
"""
function should_merge(A, B, F_prev)
  for C in F_prev
    if subseteq(A & B, C)
      return false
    end
  end
  return true
end

"""
If F contains any two sets A, B whose intersection A ∩ B is not contained in C for any C ∈ F_prev, replace A, B ∈ F with the single set A ∪ B. Repeat this operation until A ∩ B ⊆ C for some C ∈ F_prev whenever A and B are distinct members of F.

This implementation represents the sets using bits.
"""
function bitwise_superpose!(F, F_prev)
  As = copy(F)
  while length(As) !== 0
    A = pop!(As)

    for B in setdiff(F, A)
      if should_merge(A, B, F_prev)
        push!(As, A | B)
        setdiff!(F, [A, B])
        push!(F, A | B)
        break
      end
    end
  end

  return F
end

function sorted_bitwise_superpose!(F, F_prev)
  As = sort!(collect(F), by = s -> length(bits_to_set(s)))
  while length(As) !== 0
    A = popfirst!(As)

    for B in setdiff(F, A)
      if should_merge(A, B, F_prev)
        insert!(As, 1, A | B)
        setdiff!(F, [A, B])
        push!(F, A | B)
        break
      end
    end
  end

  return F
end

"""
Bitwise implementation of Knuth's approach to random matroid generation through a number of random "coarsening" steps. Supply the generate_covers and superpose methods to study the effects of different implementations of these.

n is the size of the universe.
p is a list (p_1, p_2, ...), where p_r is the number of coarsening steps to apply at rank r in the construction. The first entry of p should usually be 0, since adding closed sets of size > 1 at rank 1 is equivalent to shrinking E.
"""
function random_bitwise_kmc(generate_covers, superpose, n, p)::KnuthMatroid{Any}
  # Initialize.
  r = 1
  pr = 0
  F = [Set(0)]
  E = 2^n - 1 # The set of all elements in E.

  while true
    # Generate covers.
    push!(F, generate_covers(F[r], n))

    # Superpose.
    superpose(F[r+1], F[r])

    # Test for completion.
    if E ∈ F[r+1]
      return KnuthMatroid{Any}(n, F, [], Set(), Dict())
    end

    # Apply coarsening.
    if r <= length(p)
      pr = p[r]
    end

    while pr > 0
      # Get random closed set A in F_{r+1} and element a in E - A.
      A = rand(F[r+1])
      a = random_element(diff(E, A))

      # Replace A with A ∪ {a}.
      F[r+1] = setdiff(F[r+1], A) ∪ Set([A | a])

      # Superpose again to account for coarsening step.
      superpose(F[r+1], F[r])
      
      # Step 5: Test for completion.
      if E ∈ F[r+1]
        return KnuthMatroid{Any}(n, F, [], Set(), Dict())
      end
      
      pr -= 1
    end

    r += 1
  end
end

"""
Second implementation of random-KMC. This uses the bit-based KMC methods.
"""
function random_kmc_v2(n, p, T=UInt16)
  return random_bitwise_kmc(generate_covers_v2, bitwise_superpose!, n, p)
end

"""
Third implementation of random-KMC. This sorts the sets by size before superposing.
"""
function random_kmc_v3(n, p, T=UInt16)
  return random_bitwise_kmc(generate_covers_v2, sorted_bitwise_superpose!, n, p)
end