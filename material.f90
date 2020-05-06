module Materials

    implicit none
    integer, parameter, private :: dp = 8

    contains

        subroutine linear_elasticity(matprop, C, type)
            real(dp), dimension(:), intent(in)        :: matprop
            real(dp), dimension(:,:), intent(out)     :: C
            character(len=*), optional, intent(inout) :: type
            real(dp)                                  :: E, v

            E = matprop(1)
            v = matprop(2)

            if (.not. present(type)) then
                type = "plane_stress"
            end if

            if (type == "plane_strain") then
                C = (E / ((1+v)*(1-2*v))) * reshape([1.D0 - v, v, 0.D0, v, 1.D0 - v, 0.D0, 0.D0, 0.D0, 1.D0 - 2.D0*v], [3, 3])
            else
                C = (E / (1-v**2)) * reshape([1.D0, v, 0.D0, v, 1.D0, 0.D0, 0.D0, 0.D0, (1.D0-v)/2.D0], [3, 3])
            end if

        end subroutine linear_elasticity

end module Materials
