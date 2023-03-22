"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the ground set of elements E. Set-based.
"""
function generate_covers_v1(Fr, E)
  Set([A ∪ a for A ∈ Fr for a ∈ setdiff(E, A)])
end


"""
  If F_{r+1} contains any two sets A, B whose intersection A ∩ B is not contained in C for any C ∈ F_r, replace A, B ∈ F_{r+1} by the single set A ∪ B. Repeat this operation until A ∩ B ⊆ C for some C ∈ F_r whenever A and B are distinct members of F_{r+1}.

  F and F_old should be Family: A Set of Sets of some type.

  This is the first, Set-based implementation of this method.
"""
function superpose_v1!(F, F_old)
  for A ∈ F
    for B ∈ F
      should_merge = true
      for C ∈ F_old
        if A ∩ B ⊆ C
          should_merge = false
        end
      end



      if should_merge
        setdiff!(F, [A, B])
        push!(F, A ∪ B)
      end
    end
  end

  return F
end

"""
First implementation of Knuth's random matroid construction through random "coarsening".

n is the size of the universe.
p is a list (p_1, p_2, ...), where p_r is the number of coarsening steps to apply at rank r in the construction. The first entry of p should usually be 0, since adding closed sets of size > 1 at rank 1 is equivalent to shrinking E.

This uses the Set-based KMC methods.
"""
function random_kmc_v1(n, p, T)::KnuthMatroid{Set{Integer}}
  E = Set([i for i in range(0,n-1)])
  
  # Step 1: Initialize.
  r = 1
  F = [family([])]
  pr = 0
  
  while true
    # Step 2: Generate covers.
    push!(F, generate_covers_v1(F[r], E))
    
    # Step 4: Superpose.
    superpose_v1!(F[r+1], F[r])
    
    # Step 5: Test for completion.
    if E ∈ F[r+1]
      return KnuthMatroid{Set{Integer}}(n, F, [], Set(), Dict())
    end
    
    if r <= length(p)
      pr = p[r]
    end
    
    while pr > 0
      # Random closed set in F_{r+1} and element in E \ A.
      A = rand(F[r+1])
      a = rand(setdiff(E, A))
      
      # Replace A with A ∪ {a}.
      F[r+1] = setdiff(F[r+1], A) ∪ Set([A ∪ a])
      
      # Superpose again to account for coarsening step.
      superpose_v1!(F[r+1], F[r])
      
      
      # Step 5: Test for completion.
      if E ∈ F[r+1]
        return KnuthMatroid{Set{Integer}}(n, F, [], Set(), Dict())
      end
      
      pr -= 1
    end
    
    
    r += 1
  end
end
