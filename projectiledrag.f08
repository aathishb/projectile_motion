program projectiledrag
implicit none

real::x0, x, vx0, vx                            !x values
real::y0, y, vy0, vy                            !y values
real::t, dt                                     !time values
real::angle, anglerad                           !angle values
real::ek, ep                                    !energy values
real::v                                         !velocity

real,parameter::pi=2*asin(1.0),g = 9.79,m = 5.0,b2m = 4E-5

     print*,"Enter time step"                !initialising
     read*,dt
     x0 = 0; y0 = 0; t = 0
     v = 700.0; ek = 0.5*m*(v**2); ep = 0
     print*,"Enter the angle"
     read*,angle
     anglerad = angle*(pi/180.0)
     vx0 = v*cos(anglerad)
     vy0 = v*sin(anglerad)

     open(1, file='trajectorydrag.txt'); open(2, file='vx.txt');        open(3, file='vy.txt')
     open(4, file='ek.txt');             open(5, file='ep.txt')

     write(1,*) x0, y0;                  write(2,*) t, vx0;             write(3,*) t, vy0
     write(4,*) t, ek;                   write(5,*) t, ep

     do while (y0>=0)
             call step(v, ek, ep, x0, x, vx0, vx, y0, y, vy0, vy, t, dt)
             write(1,*) x0, y0;          write(2,*) t, vx0;             write(3,*) t, vy0
             write(4,*) t, ek;           write(5,*) t, ep
     end do

contains

subroutine step(v, ek, ep, x0, x, vx0, vx, y0, y, vy0, vy, t, dt)
implicit none
real,intent(in)::dt
real,intent(inout)::v, ek, ep, x0, x, vx0, vx, y0, y, vy0, vy, t

     v = (vx0**2 + vy0**2)**0.5             !energy increments
     ek = 0.5*m*(v**2)
     ep = m*g*y0

     x = x0 + dt*vx0                        !x incrementations
     vx = vx0 - b2m*v*vx0*dt
     x0 = x
     vx0 = vx

     y = y0 + dt*vy0                        !y incrementations
     vy = vy0 - g*dt - b2m*v*vy0*dt
     y0 = y
     vy0 = vy

     t = t + dt
end subroutine

end program
