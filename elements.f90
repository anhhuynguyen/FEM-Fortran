module Quad2D

    use MeshInfo
    implicit none

    contains

        subroutine quadrature(nip, w, ip)
            integer, intent(in)                   :: nip
            real(dp), dimension(nip), intent(out) :: w, ip

            if (nip == 1) then 
                w  = [2.0D0]
                ip = [0.0D0]
            else if (nip == 2) then
                w  = [1.0D0, 1.0D0]
                ip = [-0.577350269189626D0, 0.577350269189626D0]
            else if (nip == 3) then
                w  = [0.555555555555555D0, 0.888888888888888D0, 0.555555555555555D0]
                ip = [-0.774596669241483D0, 0.00000000000000D0, 0.774596669241483D0]
            end if

        end subroutine quadrature

        subroutine jacobian(h, x, detjac, jac)
            real(dp), dimension(dim, nodesPerElement), intent(in) :: h
            real(dp), dimension(nodesPerElement, dim), intent(in) :: x
            real(dp), intent(out)                                 :: detjac
            real(dp), dimension(dim, dim), intent(out)            :: jac

            jac = matmul(h, x)
            detjac = jac(1,1) * jac(2,2) - jac(2,1) * jac(1,2)

        end subroutine jacobian

        subroutine shape_der(s, t, x, detjac, b)
            real(dp), intent(in)                                  :: s, t
            real(dp), dimension(nodesPerElement, dim), intent(in) :: x
            real(dp), intent(out)                                 :: detjac
            real(dp), dimension(dim + 1, dofsPerElt), intent(out) :: b 
            real(dp), dimension(dim, dim)                         :: invjac
            real(dp), dimension(dim, dim)                         :: adjac, jac 
            real(dp), dimension(dim, nodesPerElement)             :: g, h

            call shape_grad(s, t, h) 
            call jacobian(h, x, detjac, jac)
            adjac = reshape([jac(2,2), -jac(2,1), -jac(1,2), jac(1,1)], [2,2])
            invjac = (1 / detjac)  * adjac
            g = matmul(invjac, h)
            b = 0.0D0
            b(1, 1:dofsPerElt:dofsPerNode) = g(1, :)
            b(2, 2:dofsPerElt:dofsPerNode) = g(2, :)
            b(3, 1:dofsPerElt:dofsPerNode) = g(2, :)
            b(3, 2:dofsPerElt:dofsPerNode) = g(1, :)
            
        end subroutine shape_der
        
        subroutine shape_func(s, t, n)
            real(dp), intent(in)                              :: s, t
            real(dp), dimension(dim, dofsPerElt), intent(out) :: n 
            
            if (eltType == "QUAD4") then
                call shape_func_Q4(s, t, n)
            else if (eltType == "QUAD8") then
                call shape_fun_Q8(s, t, n)
            end if 
            
        end subroutine shape_func 
        
        subroutine shape_grad(s, t, h)
            real(dp), intent(in)                                   :: s, t
            real(dp), dimension(dim, nodesPerElement), intent(out) :: h

            if (eltType == "QUAD4") then
                call shape_grad_Q4(s, t, h)
            else if (eltType == "QUAD8") then
                call shape_grad_Q8(s, t, h)
            end if

        end subroutine shape_grad 

        subroutine shape_func_Q4(s, t, n)
            real(dp), intent(in)                              :: s, t
            real(dp), dimension(dim, dofsPerElt), intent(out) :: n 
            real(dp), dimension(nodesPerElement)              :: eltShape 
    
            eltShape = 0.25 * [(1-s)*(1-t), (1+s)*(1-t), (1+s)*(1+t), (1-s)*(1+t)]
            n = 0.0D0
            n(1, 1:dofsPerElt:dofsPerNode) = eltShape
            n(2, 2:dofsPerElt:dofsPerNode) = eltShape
    
        end subroutine shape_func_Q4
        
        subroutine shape_grad_Q4(s, t, h)
            real(dp), intent(in)                                   :: s, t
            real(dp), dimension(dim, nodesPerElement), intent(out) :: h
            h = 0.25 * transpose(reshape([t-1, 1-t, 1+t, -1-t, s-1, -1-s, 1+s, 1-s], [4,2]))
    
        end subroutine shape_grad_Q4 

        subroutine shape_fun_Q8(s, t, n)
            real(dp), intent(in)                              :: s, t
            real(dp), dimension(dim, dofsPerElt), intent(out) :: n 
            real(dp), dimension(nodesPerElement)              :: eltShape 

            eltShape(1) = 0.25 * (1 - s) * (1 - t) * (-s-t-1) 
            eltShape(2) = 0.25 * (1 + s) * (1 - t) * (s-t-1)
            eltShape(3) = 0.25 * (1 + s) * (1 + t) * (s+t-1)
            eltShape(4) = 0.25 * (1 - s) * (1 + t) * (-s+t-1)
            eltShape(5) = 0.50 * (1 - s**2) * (1 - t)
            eltShape(6) = 0.50 * (1 + s) * (1 - t**2)
            eltShape(7) = 0.50 * (1 - s**2) * (1 + t)
            eltShape(8) = 0.50 * (1 - s) * (1 - t**2)
            n = 0.0D0
            n(1, 1:dofsPerElt:dofsPerNode) = eltShape
            n(2, 2:dofsPerElt:dofsPerNode) = eltShape

        end subroutine shape_fun_Q8

        subroutine shape_grad_Q8(s, t, h)
            real(dp), intent(in)                              :: s, t
            real(dp), dimension(dim, nodesPerElement), intent(out) :: h

            h(1,1) = 0.25 * (1 - t) * (2*s + t) 
            h(1,2) = 0.25 * (1 - t) * (2*s - t)
            h(1,3) = 0.25 * (1 + t) * (2*s + t)
            h(1,4) = 0.25 * (1 + t) * (2*s - t)
            h(1,5) = s * (t - 1)
            h(1,6) = 0.50 * (1 - t**2)
            h(1,7) = -s * (1 + t)
            h(1,8) = 0.50 * (t**2 - 1)
            h(2,1) = 0.25 * (1 - s) * (2*t + s)
            h(2,2) = 0.25 * (1 + s) * (2*t - s)
            h(2,3) = 0.25 * (1 + s) * (2*t + s)
            h(2,4) = 0.25 * (1 - s) * (2*t - s)
            h(2,5) = 0.50 * (s**2 - 1)
            h(2,6) = -t * (s + 1)
            h(2,7) = 0.5 * (1 - s**2)
            h(2,8) = t * (s - 1)  

        end subroutine shape_grad_Q8
    
end module Quad2D