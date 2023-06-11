function generate_covers_v1(Fr, E)
  Set([A ∪ a for A ∈ Fr for a ∈ setdiff(E, A)])
end


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


function random_kmc_v1(n, p, T)
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
        return (E, F)
      end
      
      pr -= 1
    end
    
    
    r += 1
  end
end
