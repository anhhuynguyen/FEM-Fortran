module MeshInfo

    implicit none

    save
    integer, parameter :: dp = 8, dofsPerNode = 2, dim = 2, thickness = 1
    integer            :: numberOfElts, numberOfNodes, nodesPerElement, dofsPerElt, numberOfDofs
    character(len=5)   :: eltType

end module MeshInfo