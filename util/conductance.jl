using LinearAlgebra
using BlockDiagonals


function gc_boundary(x; gc_neighbours=ones(2), length=1.0)
    "Conductance at boundary between blocks and poms"
    
    # Area of piece-of-metal disk
    Apom = π*0.07^2/4 - π*0.05^2/4;

    # Conductance based on volume (kc_pom*min(Ac(1),Ac(2))/l_pom);
    g_pom  = x*Apom / length;   

    # Average between piece-of-metal and neighbours
    return 1/(1/g_pom + 1/gc_neighbours[1] + 1/gc_neighbours[2])
end

function conductances(x, Ns, An, gc, length_res)
    "Construct conductance matrix of HeatedRod"

    # Initialize
    K1 = zeros(Ns[1],Ns[1]);
    K2 = zeros(Ns[2],Ns[2]);
    K3 = zeros(Ns[3],Ns[3]);

    # Base conductance
    for ii = 1:Ns[1]-1 
        K1[ii:ii+1,ii:ii+1] += [-gc[1] gc[1]; gc[1] -gc[1]];
    end
    for ii = 1:Ns[2]-1 
        K2[ii:ii+1,ii:ii+1] += [-gc[2] gc[2]; gc[2] -gc[2]];
    end
    for ii = 1:Ns[3]-1 
        K3[ii:ii+1,ii:ii+1] += [-gc[3] gc[3]; gc[3] -gc[3]];
    end

    # Add linear convection
    K1 -= diagm(h_a*An[1]*ones(Ns[1]));
    K2 -= diagm(h_a*An[2]*ones(Ns[2]));
    K3 -= diagm(h_a*An[3]*ones(Ns[3]));

    # Boundary conductances of pom resistors
    gc_b1 = gc_boundary(k_12, gc_neighbours=gc[1:2], length=length_res[1])
    gc_b2 = gc_boundary(k_23, gc_neighbours=gc[2:3], length=length_res[2])

    # Insert pom resistance into conductance matrix
    K1[Ns[1],Ns[1]] -= gc_b1;
    K2[1,1]         -= gc_b1;
    K2[Ns[2],Ns[2]] -= gc_b2;
    K3[1,1]         -= gc_b2;

    # Combine blocks
    K_total = Matrix(BlockDiagonal([K1,K2,K3]))

    # Add extra pom terms in between
    K_total[Ns[1],         Ns[1]+1      ] = gc_b1
    K_total[Ns[1]+1,       Ns[1]        ] = gc_b1
    K_total[Ns[1]+Ns[2],   Ns[1]+Ns[2]+1] = gc_b2
    K_total[Ns[1]+Ns[2]+1, Ns[1]+Ns[2]  ] = gc_b2
    
    return K_total    
end