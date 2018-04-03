#!/home/german/julia4/bin/julia
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
function read_ptraj_modes(file, modes_elements, norma::Bool=true)
    modes_file = open(file, "r")
    modes_text = readdlm(modes_file, skipstart=0, skipblanks=true,
    ignore_invalid_chars=true, comments=true, comment_char='\*')
    close(modes_file)

    nmodes = modes_text[1, 5]
    ncoords = convert(Int64, modes_elements)
    lines = ceil(Int64, ncoords/7)
    rest = convert(Int64, ncoords % 7)

    eval = Array{Float64}(nmodes);
    mode = Array{Float64}(ncoords, nmodes);
    temp1 = Array{Float64}(ncoords, 1);
    temp2 = Array{Float64}(ncoords+(7-rest));

    j = lines + 1 + 2 # 1 p/ q lea la prox linea 2 por el header

    for i = 1:nmodes
        eval[i] = modes_text[j, 2]
        temp = transpose(modes_text[(j+1):(lines+j), :])
        temp2 = reshape(temp, ncoords+(7-rest))
        for k=(rest+1):7
            pop!(temp2)
        end
    mode[:, i] = temp2
        j = j + lines + 1
    end

    if norma == true
        for i=1:nmodes
            mode[: ,i] = mode[:, i] / norm(mode[:, i])
        end
    end

    return mode, eval
end
#########
function displaceAA(in_frm, in_vec, step)
    # Preparo variables
    in_top = Topology(in_frm)
    aa = convert(Int64, count_residues(in_top))
    aa_3 = aa * 3
    natoms = convert(Int64, size(in_top))
    in_xyz = positions(in_frm)

	# Determino orden de residuos (hay q actualizar el Julia Chemfiles)
    tmp = Array{Int64}(aa)
    ids = Array{Int64}(aa)
    [ ids[i+1] = convert(Int64, id((Residue(in_top, i)))) for i = 0:aa-1 ]
    idx = sortperm(ids)
    # Determino el nro de atomos de c/ aminoácido
    [ tmp[i+1] = size(Residue(in_top, i)) for i = 0:aa-1 ]
    natom_aa = tmp[idx]

  	# Adapto el vector p/ darle la misma forma q la matriz de coordenadas
	vector = Array{Float64}
	tmp_size = size(in_vec)[1]

	if tmp_size == aa_3
		vector = transpose(reshape(in_vec, 3, aa))
	elseif tmp_size == aa
		vector = in_vec
	else
		error("Input vector with wrong dimensions: ", tmp_size, "  ", (aa_3, 1))
	end

	sum_mat = Array{Float64}(3, natoms)
	cursor = 0
   	for i = 1:aa
		rango = Array{Int64}(natom_aa[i])
    	if i == 1
			sum_mat[:, 1:natom_aa[i]] = repmat(vector[1, :], 1, natom_aa[i])
			cursor = natom_aa[i]
			continue
		end
		rango = collect(cursor+1:cursor + natom_aa[i])
        sum_mat[:, rango] = repmat(vector[1, :], 1, natom_aa[i])
		cursor += natom_aa[i]
	end

   # Listo, ahora puedo mover el pdb
   out_xyz  = in_xyz + sum_mat .* step
   return out_xyz
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
        arg_type = Int
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
    "--script"
        help = "Only write output PDB"
        action = :store_true
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
println("INPDB          ", in_pdb_filename)
println("VECTOR         ", vector_filename)
println("TOP            ", top)
println("RESOLUTION     ", resolution)
println("out_pdb_filename         ", out_pdb_filename)
println("INDEX          ", parsed_args["index"])
println("script?        ", script)

# Read PDB
const in_trj = Trajectory(in_pdb_filename)
const in_frm = read(in_trj)
const in_xyz = positions(in_frm)
const aa = convert(Int64, count_residues(in_top))
const aa_3 = aa * 3
const natoms = convert(Int64, size(in_top))

in_vec = Array{Float64}
if parsed_args["index"] != 0
# Vector de PCA Amber
    try
        in_vec = read_ptraj_modes(vector_filename,
			aa_3, true)[1][:, index]
    catch
        try
            in_vec = read_ptraj_modes(vector_filename,
				natom_xyz, true)[1][:, index]
        end
    end
else
# Vector puro
	in_vec = convert(Array{Float64}, readdlm(vector_filename))
end

# In case input vector file is note found
if in_vec == Array{Float64, 1}
    error(ArgumentError(string("\n\n", vector_filename, message_vec_not_found)))
end

# Ahora desplazo
cnt = 0
if aa_3 == length(in_vec)
# El modo es de Calpha y está ordenado en una columna
    in_vec = in_vec / norm(in_vec)
    for step in -top:resolution:top
        cnt+=1
		out_filename = string(cnt, "_", out_pdb_filename)
		out_trj = Trajectory(out_filename, 'w')
		out_frm = deepcopy(in_frm)
		out_xyz = positions(out_frm)
		out_xyz = displaceAA(in_xyz, in_vec, step);
		# Y guardo
        write(out_trj, out_frm)
    end
elseif natom_xyz == length(in_vec)
# El modo es all-atom
    in_vec = in_vec / norm(in_vec)

    for step in -top:resolution:top
        out_pdb = displaceAtoms(in_xyz, in_vec, step);
        # Y guardo
        cnt+=1
        out_filename = string(cnt, "_", out_pdb_filename)
        write(out_filename, out_pdb, PDBFile)
    end

else
# El modo no tiene el tamaño adecuado
error("PDB and input vector don't match.\nPDB has ", length(in_pdb) ,
" amino acids and ", size(coordinatesmatrix(in_pdb))[1] * 3, " atoms.\nVector has ",
length(in_vec), " elements, which should correspond to ", length(in_vec) / 3, " particles.")
end


# Finalmente, hago el script. Esto va p/ casos en los q haga 1 solo
# desplazamiento
if script == true
	f = open("script_porky.py", "w")
	load = "cmd.load(\""

	write(f, "from pymol.cgo import *\n")
    write(f, "from pymol import cmd\n\n")

	write(f, "cmd.set(\"cartoon_fancy_helices\", 1)\n")
	write(f, "cmd.set(\"cartoon_transparency\", 0.5)\n")
    write(f, "cmd.set(\"two_sided_lighting\", \"on\")\n")
    write(f, "cmd.set(\"reflect\", 0)\n")
    write(f, "cmd.set(\"ambient\", 0.5)\n")
    write(f, "cmd.set(\"ray_trace_mode\",  0)\n")
    write(f, "cmd.set('''ray_opaque_background''', '''off''')\n")

	write(f, load, in_pdb_filename,"\")\n")
	write(f, load, string(cnt, "_", out_pdb_filename),"\")\n")
	write(f, load,"modevectors.py\")\n")
	write(f, "modevectors(\"", in_pdb_filename[1:end-4], "\", \"", string(cnt, "_", out_pdb_filename)[1:end-4], "\", ")
	write(f, "outname=\"modevectors\", head=0.5, tail = 0.3, cut=0.5, headrgb = \"1.0, 1.0, 0.0\", tailrgb = \"1.0, 1.0, 0.0\") ")

	close(f)
end
