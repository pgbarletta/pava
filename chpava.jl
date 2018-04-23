#!/usr/bin/env julia
###############################################################################
# utility to displace a PDB along a vector, generating several PDBs
###############################################################################
#   code by pgbarletta
#########################################################
# Safety pig included:
#
#    _._ _..._ .-',     _.._(`))
#   '-. `     '  /-._.-'    ',/
#      )         \            '.
#     / _    _    |             \
#    |  a    a    /              |
#    \   .-.                     ;
#     '-('' ).-'       ,'       ;
#        '-;           |      .'
#           \           \    /
#           | 7  .__  _.-\   \
#           | |  |  ``/  /`  /
#          /,_|  |   /,_/   /
#             /,_/      '`-'
###############################################################################
using Chemfiles
using ArgParse
##########
# functions
##########
function read_ptraj_modes(filename, modes_elements, nmodes::Int64, norma::Bool=true)
    modes_file=open(filename, "r")
    modes_text = readdlm(modes_file, skipstart=0, skipblanks=true,
    comments=true, comment_char='\*')
    close(modes_file)

    if nmodes == 0
        nmodes = modes_text[1, 5]
    end
    ncoords = convert(Int64, modes_elements)
    lines = ceil(Int64, ncoords/7)
    rest = convert(Int64, ncoords % 7)

    evals=Array{Float64}(nmodes);
    modes = Array{Float64}(ncoords, nmodes);
    temp1=Array{Float64}(ncoords, 1);
    temp2 = Array{Float64}(ncoords+(7-rest));

    j=lines + 1 + 2 # 1 p/ q lea la prox linea 2 por el header

    for i = 1:nmodes
        evals[i] = modes_text[j, 2]
        temp = permutedims(modes_text[(j+1):(lines+j), :], [2, 1])
        temp2 = reshape(temp, ncoords+(7-rest))
        for k=(rest+1):7
            pop!(temp2)
        end
        modes[:, i] = temp2
        j = j + lines + 1
    end

    if norma == true
        for i=1:nmodes
            modes[: ,i] = modes[:, i] / norm(modes[:, i])
        end
    end

    return modes, evals
end
#########
function displaceAA(in_frm, aa, aa_3, in_vec)
    # Preparo variables
    const in_top = Topology(in_frm)
    natoms = convert(Int64, size(in_top))
    const in_xyz = positions(in_frm)

    # Determino orden de residuos (hay q actualizar el Julia Chemfiles)
    tmp = Array{Int64}(aa)
    ids = Array{Int64}(aa)
    [ ids[i+1] = convert(Int64, id((Residue(in_top, i)))) for i = 0:aa-1 ]
    const idx = sortperm(ids)
    # Determino el nro de atomos de c/ aminoácido
    [ tmp[i+1] = size(Residue(in_top, i)) for i = 0:aa-1 ]
    const natom_aa = tmp[idx]

    # Paso el vector columna de tamaño 1xaa_3 a 3xaa
    const vector = reshape(in_vec, 3, aa)
    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    sum_mat = Array{Float64}(3, natoms)
    cursor = 0
    for i = 1:aa
        if i == 1
            sum_mat[:, 1:natom_aa[i]] = repmat(vector[:, 1], 1, natom_aa[i])
            cursor = natom_aa[i]
            continue
        end
        const rango = collect(cursor+1:cursor + natom_aa[i])
        sum_mat[:, rango] = repmat(vector[:, i], 1, natom_aa[i])
        cursor += natom_aa[i]
    end

    # Listo, ahora puedo mover el pdb
    out_frm = deepcopy(in_frm)
    out_xyz = positions(out_frm)

    # Tengo q hacer esto por ahora. Hasta q arreglemos Chemfiles.
    for i = 1:size(in_xyz)[1]
        for j = 1:size(in_xyz)[2]
            out_xyz[i, j]  = in_xyz[i, j] + sum_mat[i, j]
        end
    end
    return out_frm
end
#########
function displaceAtoms(mod_pdb, vector1, multiplier)
  # Preparo variables
    pdb = copy(mod_pdb)
    struct_xyz = coordinatesmatrix(pdb)
#    new_struct_xyz = copy(struct_xyz)
    vector = Array{Float64}(1, 3)

    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    for i = 1:3:length(vector1)
        if i== 1
            vector = reshape(vector1[i:i+2], 1, 3)
            continue
        end
        vector = vcat(vector, reshape(vector1[i:i+2], 1, 3))
    end

    # Listo, ahora puedo mover el pdb
    new_struct_xyz  = struct_xyz + vector .* multiplier
    pdb = change_coordinates(pdb, new_struct_xyz);
   return pdb
end
#########
# Arg Parse settings
s = ArgParseSettings()
@add_arg_table s begin
    "--in_pdb_filename", "-p"
        help = "Input PDB"
        arg_type = String
        required = true
    "--vector_filename", "-v"
        help = "Input vector"
        arg_type = String
        required = true
    "--top", "-t"
        help = "Top/bottom multiplier"
        arg_type = Float64
        required = true
    "--resolution", "-r"
        help = "Step size between output structures. Along with --top parameter
determines the number of output structures "
        arg_type = Float64
        required = true
    "--out_pdb_filename", "-o"
        help = "Output PDBs suffix"
        arg_type = String
        required = true
    "--index", "-i"
        help = "Mode number when reading cpptraj PCA modes"
        arg_type = Int
        default = 0
end

##########
# main program
##########

# Read arguments from console
parsed_args = parse_args(ARGS, s)
args = Array{Any, 1}(0)
for (arg, val) in parsed_args
    arg = Symbol(arg)
    @eval (($arg) = ($val))
end
# Append ".pdb" to output pdb
out_pdb_filename = out_pdb_filename * ".pdb"

# Check --top argument is actually bigger than step size
if top < resolution
    throw(ArgumentError(string("\n\n --top argument ", top," is smaller/equal
than --resolution argument ", resolution)))
end

println("Input parameters:")
println("INPDB            ", in_pdb_filename)
println("VECTOR           ", vector_filename)
println("TOP              ", top)
println("RESOLUTION       ", resolution)
println("out_pdb_filename ", out_pdb_filename)
println("INDEX            ", index)

# Read PDB
const in_trj = Trajectory(in_pdb_filename)
const in_frm = read(in_trj)
const in_top = Topology(in_frm)
const aa = convert(Int64, count_residues(in_top))
const aa_3 = aa * 3
const natoms = convert(Int64, size(in_top))
const natoms_xys = natoms * 3

in_vec = Array{Float64, 1}
if index != 0
# Vector de PCA Amber
    try
        in_vec = read_ptraj_modes(vector_filename,
			aa_3, 0, true)[1][:, index]
    catch e0
        try
            in_vec = read_ptraj_modes(vector_filename,
				natom_xyz, 0, true)[1][:, index]
        catch e1
            error(e0)
            error(e1)
        end
    end
else
# Vector puro
    try
	   in_vec = convert(Array{Float64, 1}, readdlm(vector_filename))
    catch e
        error(ArgumentError(string(vector_filename,
            " could not be found.")))
    end
end

# Ahora desplazo
cnt = 0
if aa_3 == length(in_vec)
# El modo es de Calpha y está ordenado en una columna
    in_vec = in_vec / norm(in_vec)

    for step in -top:resolution:top
        cnt+=1
		const out_filename = string(cnt, "_", out_pdb_filename)
	    const out_frm = displaceAA(in_frm, aa, aa_3, in_vec .* step);
	    pos = positions(out_frm)
	    out_trj = Trajectory(out_filename, 'w')
	    write(out_trj, out_frm)
    end
elseif natom_xyz == length(in_vec)
    # El modo es all-atom
    in_vec = in_vec / norm(in_vec)

    for step in -top:resolution:top
		cnt += 1
        const out_filename = string(cnt, "_", out_pdb_filename)
		const out_frm = displaceAtoms(in_frm, aa, aa_3, in_vec .* step);
		pos = positions(out_frm)
        out_trj = Trajectory(out_filename, 'w')
        write(out_trj, out_frm)
    end
else
    # El modo no tiene el tamaño adecuado
    error("PDB and input vector don't match.\nPDB has ", length(in_pdb) ,
    " amino acids and ", size(coordinatesmatrix(in_pdb))[1] * 3, " atoms.\nVector has ",
    length(in_vec), " elements, which should correspond to ", length(in_vec) / 3, " particles.")
end
