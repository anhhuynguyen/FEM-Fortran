module Parser

    implicit none
    integer, parameter, private :: dp = 8, dim = 2

    contains

        subroutine read_input(inputFile, nodes, elements, mat, cload, boundary)

            character(len=*), intent(in)                         :: inputFile
            real(dp), dimension(:), allocatable, intent(out)     :: mat
            real(dp), dimension(:,:), allocatable, intent(out)   :: nodes, cload, boundary ! add pload for pressure
            integer, dimension(:,:), allocatable, intent(out)    :: elements
            integer                                              :: numberOfNodes, numberOfElements, nodesPerElement, &
                                                                    constrainedNodes, loadedNodes, status, msg, i
            character(len=20)                                    :: keyword, eltType, matType

            open(unit=1, file=inputFile, status='OLD', action='READ', iostat=status)

            readif: if (status == 0) then

                readloop: do 

                    read(1,*,iostat=status)  keyword
                    keyword = trim(keyword)

                    if (status /= 0) exit

                    keywordif: if (keyword == "*INFO") then
                        read(1,*) numberOfNodes, numberOfElements, constrainedNodes, loadedNodes, eltType
                        if (eltType == "QUAD4") then
                            nodesPerElement = 4
                        else if (eltType == "QUAD8") then
                            nodesPerElement = 8
                        end if
                        
                    else if (keyword == "*NODE") then
                        allocate(nodes(numberOfNodes, dim), stat=msg)
                        if (msg == 0) then
                            do i = 1, numberOfNodes
                                read(1,*) nodes(i,:)
                            end do
                        end if

                    else if (keyword == "*ELEMENT") then
                        allocate(elements(numberOfElements, nodesPerElement), stat=msg)
                        if (msg == 0) then
                            do i = 1, numberOfElements
                                read(1,*) elements(i,:)
                            end do
                        end if

                    else if (keyword == "*MATERIALS") then
                        read(1,*) matType
                        if (matType == "Elastic") then 
                            allocate(mat(2), stat=msg)
                            read(1,*) mat(:)
                        end if

                    else if (keyword == "*BOUNDARY") then
                        allocate(boundary(constrainedNodes, 5), stat=msg)
                        if (msg == 0) then
                            do i = 1, constrainedNodes
                                read(1, *) boundary(i,:)
                            end do
                        end if

                    else if (keyword == "*CLOAD") then
                        allocate(cload(loadedNodes, 3), stat=msg)
                        if (msg == 0) then
                            do i = 1, loadedNodes
                                read(1, *) cload(i,:)
                            end do
                        end if

                    end if keywordif

                end do readloop

            else 
                print *, "Failed to open file"

            end if readif

            close(1)

        end subroutine read_input

        subroutine check_input(nodes, elements, mat, disp_bc)
            integer, dimension(:,:), allocatable, intent(in)  :: elements
            real(dp), dimension(:,:), allocatable, intent(in) :: nodes, disp_bc
            real(dp), dimension(:), allocatable, intent(in)   :: mat

            if (.not. allocated(elements) .or. .not. allocated(nodes) .or. &
                .not. allocated(mat) .or. .not. allocated(disp_bc)) then 
                    print*, "Something not allocated"
                call exit() 
            end if

        end subroutine check_input

        subroutine write2file(nodes, u, sigma, eps, inputFile, outputFile)
            real(dp), dimension(:,:), intent(in) :: u, sigma, eps, nodes
            character(len=*), intent(in)         :: outputFile, inputFile
            integer                              :: status, nip, i, j, k
            character(len=30)                    :: date 

            if (size(sigma,2) == 3) then
                nip = 1
            else if (size(sigma,2) == 12) then
                nip = 2
            else if (size(sigma,2) == 27) then
                nip = 3
            end if
            
            open(unit=2, file=outputFile, status="REPLACE", action="WRITE", iostat=status)
            
            100 format (A15, 1X, A10, 1X, A15, 1X, A15, 1X, A15)
            101 format (I15, 1X, I10, 1X, ES15.5, 1X, ES15.5, 1X, ES15.5)
            102 format (A15, 1X, A15, 1X, A15, 1X, A15, 1X, A15, 1X, A15)
            103 format (I15, 1X, ES15.5, 1X, ES15.5, 1X, ES15.5, 1X, F15.5, 1X, F15.5)

            call fdate(date)

            writeif: if (status == 0) then

                write (2,*) repeat("*", 94) 
                write (2,*) "Field Output Report, written at ", date
                write (2,*) 
                write (2,*) "Input File: ", inputFile
                write (2,*) 


                write (2,*) repeat('-', 37), "Nodal Displacements", repeat('-', 38)
                write (2,102) "Node Label", "Magnitude", "U1", "U2", "X", "Y" 
                write (2,*) repeat("-", 94)
                
                disploop: do i = 1, size(u, 1), size(nodes, 2)
                    write (2,103) ceiling(i/2.), sqrt(sum(u(i:i+1,1)**2)), u(i,1), u(i+1,1), nodes(ceiling(i/2.),:)
                end do disploop
                
                write (2,*) repeat("-", 94)
                
                write (2,*) 
                write (2,*) repeat('-', 22), "Stresses at integration points", repeat('-', 21)
                write (2,100) "Element Label", "IntPt", "Sigma_11", "Sigma_22", "Sigma_12"
                write (2,*) repeat("-", 73)

                
                stressloop: do i = 1, size(sigma, 1)
                k = 1
                do j = 1, nip**2
                    write (2,101) i, j, sigma(i, k:k+2)
                    k = k + 3
                end do 
                end do stressloop
                
                write (2,*) repeat("-", 73)
            
                write (2,*) 
                write (2,*) repeat('-', 22), "Strains at integration points", repeat('-', 22)
                write (2,100) "Element Label", "IntPt", "Eps_11", "Eps_22", "Eps_12"
                write (2,*) repeat("-", 73)
                
                strainloop: do i = 1, size(eps, 1)
                k = 1
                do j = 1, nip**2
                    write (2,101) i, j, eps(i, k:k+2)
                    k = k + 3
                end do 
                end do strainloop
                
                write (2,*) repeat("-", 73)
            
            else 
                write (*,*) "Fail to write solutions"

            end if writeif

            close(2)

        end subroutine write2file

end module Parser




