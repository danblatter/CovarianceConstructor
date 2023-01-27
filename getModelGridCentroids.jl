
module getModelGridCentroids  # this is a module so variables defined outside the "for" loop
# become global in scope and are thus available inside the for
# loop. Otherwise Julia complains A LOT!

# imports:
import MPI, Mare2dem
using DelimitedFiles

# MARE2DEM parameter files:
filenameroot = "KI_60_justReservoir2.0" # MARE2DEM model files to load
bquiet       = true       # set to true to turn off all MARE2DEM print statements
outfilename = "KI_60centroids_justReservoir2.txt"
outfilename2 = "KI_60rho_justReservoir2.txt"

# MPI stuff
# Initialize MPI and get communicator, rank and number of MPI processes in communicator:
MPI.Init()

comm  = MPI.COMM_WORLD           # the default MPI communicator, which has all processes invoked with mpirun
mpirank  = MPI.Comm_rank(comm)     # the id (index) of this process in communicator comm
nproc = MPI.Comm_size(comm)     # the total number of processes in the communicator group comm

if mpirank == 0
    M2d = Mare2dem.load(filenameroot,bquiet)
    # print field names:
    println("field names in returned struct:")
    for fname in fieldnames(typeof(M2d))
        println("$fname ")
    end
    writedlm(outfilename,M2d.centroids)
    writedlm(outfilename2,M2d.log10rhofree)
end

MPI.Barrier(comm)
MPI.Finalize()

end