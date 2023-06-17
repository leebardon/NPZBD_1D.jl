using Parameters

H = 400          
dz = 10                       
ngrid = Int(H/dz)             
zc = [dz/2 : dz : H - dz/2;]   
zf = [0 : dz : H;]      

nd = 3
ws = zeros(nd)                  
ws[1] = 10 
w = zeros(ngrid + 1)            
wd = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd)
wd[1,:] .= 0                   
wd[end,:] .= 0  

mlz=25                      
kappazmin=1e-4            
kappazmax=1e-2              
kappa_z=(kappazmax .* exp.(-zf/mlz) .+ kappazmin .+ kappazmax .* exp.((zf .- H) / 100.)) .* 3600 .* 24 
kappa_z[1]=0
kappa_z[end]=0

euz = 25                   
light_top = 700             
light = light_top .* exp.( -zc ./ euz)   

temp = 12 .*exp.(-zc./ 150) .+ 12 .*exp.(-zc ./ 500) .+ 2
temp_coeff_arr = 0.8
temp_ae_arr = -4000
temp_ref_arr = 293.15   
t_kel = 273.15
temp_fun = temp_coeff_arr .* exp.(temp_ae_arr .*(1 ./ (temp .+ t_kel) .- 1 ./ temp_ref_arr))


@with_kw struct TestPrms
    tt::Int64  = 2
    dt::Float64  = 0.01                  
    nrec::Int64  = 100 
    H::Int64  = 400                 
    dz::Int64  = 10              
    np::Int64    = 1                 
    nb::Int64    = 6                
    nz::Int64    = 2                  
    nn::Int64    = 1                 
    nd::Int64    = nd                 
    pIC::Array{Float64,2} = ones(Float64, 40, 1) * 0.1     
    bIC::Array{Float64,2} = ones(Float64, 40, 6) * 0.01         
    zIC::Array{Float64,2} = ones(Float64, 40, 2) * 0.01         
    nIC::Array{Float64,2} = ones(Float64, 40, 1) * 5.0        
    dIC::Array{Float64,2} = ones(Float64, 40, 3) * 0.1      
    oIC::Array{Float64,2} = ones(Float64, 40, 1) * 100.0  
    umax_p::Array{Float64,1} = [1.0]     
    K_n::Array{Float64,1}  = [0.1]       
    m_lp::Array{Float64,1} = [0.1]      
    m_qp::Array{Float64,1} = [1.0]    
    CM::Array{Bool,2}  = [1 0 0 1 0 0; 0 1 0 0 0 1; 0 0 1 0 1 0]          
    y_ij::Array{Float64,2} = [0.3 0.0 0.0 0.3 0.0 0.0; 0.0 0.3 0.0 0.0 0.0 0.3; 0.0 0.0 0.3 0.0 0.3 0.0]       
    vmax_i::Array{Float64,1} = [0.010000000000000002, 0.31622776601683794, 10.0]    
    vmax_ij::Array{Float64,2} = [0.0180815 0.0 0.0 0.0 0.0 0.0; 0.0 0.446873 0.0 0.0 0.161262 0.0; 0.0 0.0 9.40267 12.8549 0.0 16.5146]   
    Km_i::Array{Float64,1}  = [0.0010000000000000002, 0.03162277660168379, 1.0]       
    Km_ij::Array{Float64,2} = [0.000335257 NaN NaN 0.000622884 NaN NaN; NaN 4.27197 NaN NaN NaN 0.160194; NaN NaN 2.51638 NaN 2.00226 NaN]      
    m_lb::Array{Float64,1}  = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]       
    m_qb::Array{Float64,1}  = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]       
    g_max::Array{Float64,1} = [1.0, 1.0]      
    K_g::Array{Float64,1} = [0.2, 0.2]       
    γ::Array{Float64,1}  = [0.3, 0.3]          
    m_lz::Array{Float64,1} = [0.01, 0.01]       
    m_qz::Array{Float64,1} = [1.0, 1.0]       
    fsave::String  = "test_out_1D"                  
    GrM::Array{Bool,2}  = [0 1 1 1 1 1 1; 1 0 1 0 1 1 1]           
    pen::Array{Float64,1} = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]         
    SW::Array{Float64,1} = [1.0, 1.0, 1.0]      
    prob_generate_d::Array{Float64, 1} = [0.3, 0.6, 0.1]
    kappa_z::Array{Float64, 1} = kappa_z
    wd::Array{Float64, 2} = wd
    light::Array{Float64, 1} = light
    temp_fun::Array{Float64, 1} = temp_fun
    K_I::Array{Float64, 1} = [10.0]
    ngrid::Int64 = 40
    e_o::Float64 = 9.375
    yo_ij::Array{Float64, 2} = [3.0 0.0 0.0 3.0 0.0 0.0; 0.0 3.0 0.0 0.0 0.0 3.0; 0.0 0.0 3.0 0.0 3.0 0.0]  
    koverh::Float64 = 0.002592
    o2_sat::Float64 = 212.1
    ml_boxes::Int64 = 10
    t_o2relax::Float64 = 0.01
    o2_deep::Float64 = 200.0     
end
