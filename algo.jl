using Oscar;
# Temporary: Using a fixed example to ensure the code works, will generalise to n parameters soon
# We write Ka for K[a]
Ka, (a1,a2,a3,a4,a5,a6,a7,a8) = polynomial_ring(QQ,["a1","a2","a3","a4","a5","a6","a7","a8"])
S,(x,y,z1,z2,z3,z4,z5,z6,z7,z8)=polynomial_ring(Ka,["x","y","z1","z2","z3","z4","z5","z6","z7","z8"])
f=[a1+a2*x^2+a3*y^5+a4*(x^2+y^2+x^2*y^3),a5+a6*x+a7*y+a8*(x^2+y^2)]

# using findfirst each time is too much typing

function get_position(element,list)
    n=findfirst(x->x==element,list)
    return n
end

function get_polynomial_for_coeff(del,pol)
    phi=hom(S,S,c->evaluate(c,[Int(a==del) for a in gens(Ka)]),gens(S))
    return phi(pol)
end

# a list of terms associated to each parameter (a1,...) per polynomial
hf=[[] for g in f]

for g in f
    for j in gens(Ka)
        push!(hf[get_position(g,f)],get_polynomial_for_coeff(j,g))
    end
end

# strip repeated monomials and 0
hflat=vcat(hf...)
hred=unique(hflat)
hred=filter!(!=(0),hred)
# we only replace non-monomial terms
hred=filter!(x->~is_monomial(x),hred)

mod_system=[]
for g in f
    mod_f=0
    for j in hf[get_position(g,f)]
        if j in hred
            # replace non monomial term with variable, add j_i -z_i to modified system
            mod_f+=gens(Ka)[get_position(j,hf[get_position(g,f)])]*gens(S)[get_position(j,hred)+2]
            push!(mod_system,j-gens(S)[get_position(j,hred)+2])
        else
            mod_f+=gens(Ka)[get_position(j,hf[get_position(g,f)])]*j
        end
    end
    push!(mod_system,mod_f)
end
#strip repeats
mod_system=unique(mod_system)

gen_coeff=hom(S,S,c->evaluate(c,[(Int(rand(Int16))) for a in gens(Ka)]),gens(S))
# TODO: map these to t^(c) so they can be sent directly to gfan
for g in mod_system
    println(map_coefficients(c->c, gen_coeff(g)))
end

    
