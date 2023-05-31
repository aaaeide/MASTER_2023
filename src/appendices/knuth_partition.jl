"""
Knuth's 1973 Matroid Partitioning algorithm for partitioning a set into subsets 
independent in various given matroids.

Knuth's description: Given k matroids Ms = [M1, ..., Mk] on the same ground set 
E, the algorithm finds a k-partitioning [S1, ..., Sk] of the elements of E such 
that Sj is independent in matroid Mj, and nj <= |Sj| <= nj', for given limits 
nj and nj'. 

This implementation drops the upper limit nj' for each element j (implicitly it 
is infinity for all matroids). Supply nj in array lims (lims[j] = nj).
"""
function knuth_partition(Ms, lims=nothing)
  n = Ms[1].n; k = length(Ms)
  S0 = Set(1:n) # The unallocated items.
  S = [Set() for _ in 1:k] # The partition-to-be.
  color = Dict(x=>0 for x in 1:n) # color[x] = j iff x ∈ S[j].
  for y in 1:k color[-y] = y end # -y is the 'standard' element of color y.
  succ = [0 for _ in 1:n]

  lims = lims === nothing ? [0 for _ in 1:k] : lims

  function augment(r)
    for x in 1:n succ[x] = 0 end
    
    A = Set(1:n)
    B = r > 0 ? Set(-r) : Set(-j for j in 1:k)
    
    while B != Set()
      C = Set()
      for y ∈ B for x ∈ A
        j = color[y]

        if is_indep(Ms[j], x ∪ setdiff(S[j], y))
          succ[x] = y
          A = setdiff(A, x)
          C = C ∪ x
          if color[x] == 0 repaint(x); return end
        end
      end end
      B = C
    end

    println("$A violates the condition of Theorem 3")
  end

  function repaint(x)
    while x ∈ 1:n
      y = succ[x]
      j = color[x]
      
      if j == 0 setdiff!(S0, x) else setdiff!(S[j], x) end
      
      j = color[y]
      S[j] = S[j] ∪ x
      color[x] = j
      x = y
    end
  end

  # Ensure every part gets at least its lower limit.
  for j in 1:k for _ in 1:lims[j]
    augment(j)
  end end

  # Allocate the rest.
  while S0 != Set()
    augment(0)
  end

  return S
end