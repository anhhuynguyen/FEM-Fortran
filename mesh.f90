module MeshInfo

    implicit none

    save
    integer, parameter :: dp = 8, dofsPerNode = 2, dim = 2, thickness = 1
    integer            :: numberOfElts, numberOfNodes, nodesPerElement, dofsPerElt, numberOfDofs
    character(len=5)   :: eltType

    interface printMatrix
        module procedure printMatrix_i
        module procedure printMatrix_d
        module procedure printVector
    end interface

    contains

    subroutine print_matrix(arr)
        real(dp), intent(in), dimension(:,:) :: arr
        integer                              :: i, j

        do i = 1, size(arr, 1)
            write (*,100) (arr(i, j), j = 1, size(arr, 2))
            100 format (30 E15.4)
        end do

    end subroutine print_matrix

    subroutine printMatrix_i(array)
        implicit none
        integer, intent(in) :: array(:,:)
        integer :: i
        
        do i = 1, size(array, 1)
            print*, array(i,:)
        end do

    end subroutine printMatrix_i
    
    subroutine printMatrix_d(array)
        implicit none
        real(dp), intent(in) :: array(:,:)
        integer :: i
        
        do i = 1, size(array, 1)
            print*, array(i,:)
        end do

    end subroutine printMatrix_d
    
    subroutine printVector(array)
        implicit none
        real(dp), intent(in) :: array(:)
        integer :: i
        
        do i = 1, size(array, 1)
            print*, array
        end do

    end subroutine printVector

end module MeshInfo