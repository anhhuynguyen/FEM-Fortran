module Solver

    use Parser
    use MeshInfo 
    use Quad2D
    use Materials

    implicit none
    
    contains 

        subroutine initialize(nodes, elements, k, f, ipiv)
            real(dp), dimension(:,:), intent(in)               :: nodes
            integer, dimension(:,:), intent(in)                :: elements
            real(dp), dimension(:,:), allocatable, intent(out) :: k, f
            integer, dimension(:), allocatable, intent(out)    :: ipiv
            integer                                            :: k_msg, f_msg

            numberOfElts = size(elements, 1)
            nodesPerElement = size(elements, 2)
            numberOfNodes = size(nodes, 1)
            dofsPerElt = nodesPerElement * dofsPerNode 
            numberOfDofs = numberOfNodes * dofsPerNode

            if (nodesPerElement == 4) then
                eltType = "QUAD4"
            else if (nodesPerElement == 8) then
                eltType = "QUAD8"
            end if

            allocate(k(numberOfDofs, numberOfDofs), STAT=k_msg)
            allocate(f(numberOfDofs, 1), STAT=f_msg)
            allocate(ipiv(numberOfDofs))

            k = 0.
            f = 0.

            if (k_msg == 0 .and. f_msg == 0) then
                print *, "Succeed to alloacate memory"
            else
                print *, "Fail to allocate memory. Stop solving"
            end if

        end subroutine initialize

        subroutine local_to_global_dofs(element, dofs)
            integer, dimension(:), intent(in)  :: element
            integer, dimension(:), intent(out) :: dofs
            integer                            :: j
            dofs = [( dofsPerNode*element(j) - 1, dofsPerNode*element(j), j = 1, size(element) )]
            
        end subroutine local_to_global_dofs

        subroutine stiffness(c, x, nip, k)
            real(dp), dimension(:,:), intent(in)              :: c, x
            integer , intent(in)                              :: nip
            real(dp), dimension(:,:), intent(out)             :: k
            real(dp), dimension(nip)                          :: w, ip
            real(dp), dimension(dofsPerNode + 1, dofsPerElt)  :: b
            real(dp)                                          :: detjac 
            integer                                           :: i, j 

            call quadrature(nip, w, ip)

            k = 0.
            do i = 1, nip
                do j = 1, nip
                    call shape_der(ip(i), ip(j), x, detjac, b)
                    k = k + matmul(matmul(transpose(b), c), b) * detjac * w(i) * w(j) * thickness
                end do
            end do

        end subroutine stiffness

        subroutine stiffness_assembler(elements, nodes, mat, k)
            real(dp), dimension(:,:), intent(in)        :: nodes
            integer, dimension(:,:), intent(in)         :: elements
            real(dp), dimension(:,:), intent(in)        :: mat
            real(dp), dimension(:,:), intent(out)       :: k
            real(dp), dimension(dofsPerElt, dofsPerElt) :: kel
            integer, dimension(nodesPerElement)         :: element
            integer, dimension(dofsPerElt)              :: dofs
            integer :: i, m, n, nip

            if (eltType == "QUAD4") then
                nip = 2
            else if (eltType == "QUAD8") then
                nip = 3
            end if 
            
            do i = 1, numberOfElts
                element = elements(i,:)
                call stiffness(mat, nodes(element,:), nip, kel)
                call local_to_global_dofs(element, dofs)
                do m = 1, size(dofs)
                    do n = 1, size(dofs)
                        k(dofs(m), dofs(n)) = k(dofs(m), dofs(n)) + kel(m, n)
                    end do
                end do
            end do

        end subroutine stiffness_assembler

        subroutine force_assembler(force, f)
            real(dp), dimension(:,:), allocatable, intent(in)  :: force
            real(dp), dimension(:,:), intent(out)              :: f
            integer                                            :: i, j

            ! Think about vectorization
            ! point loads => call apply_poit_load()
            if (allocated(force)) then 
                do i = 1, size(force, 1)
                    j = int(force(i, 1))
                    f(dofsPerNode*j - 1, 1) = f(dofsPerNode*j - 1, 1) + force(i, 2)
                    f(dofsPerNode*j, 1)   = f(dofsPerNode*j, 1) + force(i, 3)
                end do
            end if
            ! surface loads
            ! call apply_pressure()

        end subroutine force_assembler

        subroutine disp_imposer(disp, k, f)
            real(dp), dimension(:,:), intent(in)    :: disp
            real(dp), dimension(:,:), intent(inout) :: k, f
            integer                                 :: i, j, constraint_x, constraint_y
            real(dp), parameter                     :: penaltyFactor = 1.0e20
            real(dp)                                :: dx, dy

            do i = 1, size(disp, 1)
                j = int(disp(i, 1))
                dx = disp(i, 2)
                dy = disp(i, 3)
                constraint_x = int(disp(i, 4))
                constraint_y = int(disp(i, 5))

                if (constraint_x == 1) then
                    k(dofsPerNode*j - 1, dofsPerNode*j - 1) = k(dofsPerNode*j - 1, dofsPerNode*j - 1) * penaltyFactor
                    f(dofsPerNode*j - 1, 1) = k(dofsPerNode*j - 1, 2*j - 1) * dx
                end if
                
                if (constraint_y == 1) then
                    k(dofsPerNode*j, dofsPerNode*j) = k(dofsPerNode*j, dofsPerNode*j) * penaltyFactor
                    f(dofsPerNode*j, 1) = k(dofsPerNode*j, dofsPerNode*j) * dy
                end if
            end do

        end subroutine disp_imposer

        subroutine get_stress(elements, nodes, nip, u, c, epsilon, sigma)
            real(dp), dimension(:,:), intent(in)                 :: u, c, nodes
            integer, dimension(:,:), intent(in)                  :: elements
            integer, intent(in)                                  :: nip
            real(dp), dimension(:,:), allocatable, intent(inout) :: sigma, epsilon
            real(dp), dimension(dofsPerNode + 1, dofsPerElt)     :: b
            real(dp), dimension(nip)                             :: w, ip
            real(dp)                                             :: detjac 
            integer, dimension(nodesPerElement)                  :: element
            integer, dimension(dofsPerElt)                       :: dofs
            integer                                              :: el, i, j, k

            allocate(sigma(numberOfElts, 3*nip**2))
            allocate(epsilon(numberOfElts, 3*nip**2))
            
            call quadrature(nip, w, ip)

            elementloop: do el = 1, numberOfElts
                element = elements(el,:)
                k = 1
                do i = 1, nip
                    do j = 1, nip
                        call shape_der(ip(j), ip(i), nodes(element,:), detjac, b)
                        call local_to_global_dofs(element, dofs)
                        call compute_strain(b, u(dofs,1), epsilon(el, k:k+2))
                        call compute_stress(c, epsilon(el, k:k+2), sigma(el, k:k+2))
                        k = k + 3
                    end do
                end do
            end do elementloop

        end subroutine get_stress

        subroutine compute_strain(b, u, eps)
            real(dp), dimension(:,:), intent(in) :: b
            real(dp), dimension(:), intent(in)   :: u
            real(dp), dimension(:), intent(out)  :: eps

            eps = matmul(b, u)
        
        end subroutine compute_strain

        subroutine compute_stress(c, eps, sigma)
            real(dp), dimension(:,:), intent(in) :: c
            real(dp), dimension(:), intent(in)   :: eps
            real(dp), dimension(:), intent(out)  :: sigma

            sigma = matmul(c, eps) 

        end subroutine compute_stress

        subroutine stiffness_destructor(k)
            real(dp), dimension(:,:), allocatable, intent(inout) :: k
            integer                                              :: msg

            if (allocated(k)) then
                deallocate(k, stat=msg)
                if (msg /= 0 ) then
                    write(*,*) "Fail to deallocate the global stiffness matrix"
                end if
            end if 

        end subroutine stiffness_destructor

        subroutine solve(inputfile, outputfile)
            character(len=*), intent(in)          :: inputfile, outputfile
            real(dp), dimension(:,:), allocatable :: nodes, force_bc, disp_bc
            integer, dimension(:,:), allocatable  :: elements
            real(dp), dimension(:), allocatable :: mat
            real(dp), dimension(dim+1,dim+1)      :: c
            real(dp), dimension(:,:), allocatable :: k, f, stress, strain
            integer, dimension(:), allocatable    :: ipiv
            integer                               :: info

            call read_input(inputfile, nodes, elements, mat, force_bc, disp_bc)
            call check_input(nodes, elements, mat, disp_bc)
            call initialize(nodes, elements, k, f, ipiv)
            call linear_elasticity(mat, c)
            call stiffness_assembler(elements, nodes, c, k)
            call force_assembler(force_bc, f)
            call disp_imposer(disp_bc, k, f)
            call dgesv(numberOfDofs, 1, k, numberOfDofs, ipiv, f, numberOfDofs, info)
            call stiffness_destructor(k)
            call get_stress(elements, nodes, 1, f, c, strain, stress)
            call write2file(nodes, f, stress, strain, inputfile, outputfile)

        end subroutine solve

end module Solver
