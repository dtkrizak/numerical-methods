!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module importantEqs saves the important variables/equations to be used by both the
! program and the subroutine
Module importantEqs

!The functions are added into the module for use in the subroutines
Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function to be solved, the equation for the valve for the tank with
!
Function func(t, F1, tau)
Double Precision, intent (in):: t, F1
Double Precision:: func, F1spec, tau
!tau = 2.00
!F1spec is a function of time F1(t=0) = 5kg/s]
If (t>=10.0) then
F1spec = 4.00
else
F1spec = 5.00 - 0.10*t
End if

!Sets the value of the function
func = 1/tau*(F1spec - F1)
End Function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function2 to be solved, the equation for the energy of the system for a tank with two
! inlets with flow of F1 (kg/s) and F2 (kg/s) at temperatures of T1 (˚C) and T2 (˚C)
! and one outlet flow (F1+F2) at temperature Temp (˚C) and an accumulation term M (kg)
Function func2(t, F1, Temp)
Double Precision, intent (in):: t, F1, Temp
Double Precision:: func2, M, F2, T1, T2
M = 100.00
T1 = 25.00
T2 = 75.00
F2 = 5.00
!Sets the value of the function
func2 = 1/M*(F1*T1 + F2*T2 - (F1+F2)*Temp)
End Function

End Module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program Proj4 utilizes different subroutines to solve for the equation for the
! temperature of the outlet stream of a tank and produces graphs of the results
Program Proj4_2
Use importantEqs
! t0 	: the starting value of time (seconds)
! tf 	: the final value of time (seconds)
! Temp0 : the starting value of temperature (˚C)
! steps : the number of steps to divide the interval
! intX	: array for x values from the Integration Method
! intY	: array for y values from the Integration Method
! rkX	: array for x values from the Runge-Kutta Method
! rkY	: array for y values from the Runge-Kutta Method
Double Precision:: t0 = 0.0, tf = 75.0, flow0 = 5.0, Temp0 = 50.0
Integer, parameter::steps = 50
Double Precision:: intX2(0:steps), intY2(0:steps), intZ2(0:steps), rkX2(0:steps), rkY2(0:steps), rkZ2(0:steps)
Double Precision:: intX4(0:steps), intY4(0:steps), intZ4(0:steps), rkX4(0:steps), rkY4(0:steps), rkZ4(0:steps)
Double Precision:: intX6(0:steps), intY6(0:steps), intZ6(0:steps), rkX6(0:steps), rkY6(0:steps), rkZ6(0:steps)
Double Precision::tau2 = 2.00, tau4 = 4.00, tau6 = 6.00
Character(len = 1)::name2 = '2', name4 = '4', name6 = '6'

!Uses the variables from the module to solve func for both methods
Call IntegrationMethod(t0, tf, flow0, Temp0, steps, intX2, intY2, intZ2, tau2, name2)
Call RungeKuttaMethod(t0, tf, flow0, Temp0, steps, rkX2, rkY2, rkZ2, tau2, name2)
Call Combined_Data(steps, intX2, intZ2, rkZ2, name2)
Call Combined_Data(steps, intX2, intY2, rkY2, name2//'F')


!Uses the variables from the module to solve func for both methods
Call IntegrationMethod(t0, tf, flow0, Temp0, steps, intX4, intY4, intZ4, tau4, name4)
Call RungeKuttaMethod(t0, tf, flow0, Temp0, steps, rkX4, rkY4, rkZ4, tau4, name4)
Call Combined_Data(steps, intX4, intZ4, rkZ4, name4)
Call Combined_Data(steps, intX4, intY4, rkY4, name4//'F')

!Uses the variables from the module to solve func for both methods
Call IntegrationMethod(t0, tf, flow0, Temp0, steps, intX6, intY6, intZ6, tau6, name6)
Call RungeKuttaMethod(t0, tf, flow0, Temp0, steps, rkX6, rkY6, rkZ6, tau6, name6)
Call Combined_Data(steps, intX6, intZ6, rkZ6, name6)
Call Combined_Data(steps, intX6, intY6, rkY6, name6//'F')

Call run_gnuplot('/Users/DanielK/Desktop/Fortran/IntRKMethod_Part2All.txt')
Call run_gnuplot('/Users/DanielK/Desktop/Fortran/IntRKMethod_Part2AllF.txt')


End Program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine IntegrationMethod utilizes the numerical integration method to calculate
! the area under the curve for the given function in the module
Subroutine IntegrationMethod(startX, finalX, startY, startZ, numSteps, xint, yint, zint, tauin, namein)
Use importantEqs
Double Precision, intent(in):: startX, finalX, startY, startZ, tauin
Integer, intent(in)::numSteps
Character(len = 1), intent(in)::namein
Double Precision:: stepSize, k1, k2, k3, k4, F1calc
Double Precision, intent(out):: xint(0:numSteps), yint(0:numSteps), zint(0:numSteps)
Integer::err

!Calculates the steps size and sets the initial values based on the input data
stepSize = (finalX - startX)/numSteps
xint(0) = startX
yint(0) = startY
zint(0) = startZ

!Do loop runs the interative calculation of the Integration Method
Do i = 0, numSteps - 1
	xint(i+1) = xint(i) + stepSize
!	Print *, 'xint', xint(i), 'yint', yint(i)
	yint(i+1) = yint(i) + stepSize*func(xint(i), yint(i), tauin)
!	Print*, '		stepSize*func', stepSize*func(xint(i), yint(i))
	zint(i+1) = zint(i) + stepSize*func2(xint(i), yint(i), zint(i))
End Do

err = 0

!The write_xy_data subroutine writes the data in the arrays to the specified data file
Call write_xy_data('/Users/DanielK/Desktop/Fortran/IntMethod_'//namein//'.dat', numSteps, xint, zint, err)
Call write_xy_data('/Users/DanielK/Desktop/Fortran/IntMethod_'//namein//'F.dat', numSteps, xint, yint, err)

!The write_xy_plot subroutine can write the script to produce a graph of the data
!Call write_xy_plot('/Users/DanielK/Desktop/Fortran/IntMethod.txt', 'IntMethod.dat', err)

!The run_gnuplot subroutine uses the given script to create a graph of the data
!Call run_gnuplot('/Users/DanielK/Desktop/Fortran/IntMethod_2.txt')

if(err /= 0) then
Print*, 'Fatal Error: Unable to plot IntMethod_'//namein//'.dat'
End if

End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine utilizes RungeKuttaMethod utilizes the Runge-Kutta Method to calculate the
! area under the curve for the given function in the module, more detail is given at:
! https://math.stackexchange.com/questions/721076/help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od
Subroutine RungeKuttaMethod(startX, finalX, startY, startZ, numSteps, xrk, yrk, zrk, tauin, namein)
Use importantEqs
Double Precision, intent(in):: startX, finalX, startY, startZ, tauin
Integer, intent(in)::numSteps
Character(len = 1), intent(in)::namein
Double Precision:: stepSize, k1, k2, k3, k4, k1z, k2z, k3z, k4z
Double Precision, intent(out):: xrk(0:numSteps), yrk(0:numSteps), zrk(0:numSteps)
Integer::err

!Calculates the steps size and sets the initial values based on the input data
stepSize = (finalX - startX)/numSteps
xrk(0) = startX
yrk(0) = startY
zrk(0) = startZ

!Do loop runs the interative calculation of the RK Method
Do i = 0, numSteps - 1
	xrk(i+1) = xrk(i) + stepSize
!	Print *, 'xrk', xrk(i), 'yrk', yrk(i)
	k1 = stepSize*func(xrk(i), yrk(i), tauin)
	k2 = stepSize*func((xrk(i)+stepSize/2), (yrk(i)+k1/2), tauin)
	k3 = stepSize*func((xrk(i)+stepSize/2), (yrk(i)+k2/2), tauin)
	k4 = stepSize*func((xrk(i)+stepSize), (yrk(i)+k3), tauin)
!	Print *, '		k1', k1, '	k2', k2, '	k3', k3, '	k4', k4
	yrk(i+1) = yrk(i) + (k1 + 2.00*k2 + 2.00*k3 + k4)/6.00
!	Print*, '		1/6(...)', (k1 + 2.00*k2 + 2.00*k3 + k4)/6.00
	k1z = stepSize*func2(xrk(i), yrk(i), zrk(i))
	k2z = stepSize*func2((xrk(i)+stepSize/2), (yrk(i)+k1/2), (zrk(i)+k1z/2))
	k3z = stepSize*func2((xrk(i)+stepSize/2), (yrk(i)+k2/2), (zrk(i)+k2z/2))
	k4z = stepSize*func2((xrk(i)+stepSize), (yrk(i)+k3), (zrk(i)+k3z))
!	Print *, '		k1', k1, '	k2', k2, '	k3', k3, '	k4', k4
    zrk(i+1) = zrk(i) + (k1z + 2.00*k2z + 2.00*k3z + k4z)/6.00
!	Print*, '		1/6(...)', (k1 + 2.00*k2 + 2.00*k3 + k4)/6.00

End Do

err = 0

!The write_xy_data subroutine writes the data in the arrays to the specified data file
Call write_xy_data('/Users/DanielK/Desktop/Fortran/RKMethod_'//namein//'.dat', numSteps, xrk, zrk, err)
Call write_xy_data('/Users/DanielK/Desktop/Fortran/RKMethod_'//namein//'F.dat', numSteps, xrk, yrk, err)

!The write_xy_plot subroutine can write the script to produce a graph of the data
!Call write_xy_plot('/Users/DanielK/Desktop/Fortran/RKMethod.txt', 'RKMethod.dat', err)

!The run_gnuplot subroutine uses the given script to create a graph of the data
!Call run_gnuplot('/Users/DanielK/Desktop/Fortran/RKMethod_2.txt')

if(err /= 0) then
	Print*, '/Users/DanielK/Desktop/Fortran/RKMethod_'//namein//'.dat'
End if

End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine IntegrationMethod utilizes the numerical integration method to calculate
! the area under the curve for the given function in the module
Subroutine Combined_Data(numSteps, xint, yint, yrk, namein)
Use importantEqs
Integer, intent(in)::numSteps
Character(len = *), intent(in)::namein
Double Precision, intent(in):: xint(0:numSteps), yint(0:numSteps), yrk(0:numSteps)
Double Precision::allpts(0:2, 0:numSteps)
Integer::err

!Do loops inputs all the data from both solutions into a single array
Do i = 0, numSteps
	allpts(i,0) = xint(i)
	allpts(i,1) = yint(i)
	allpts(i,2) = yrk(i)
!	Print*, 'allpts(0,i)', allpts(0,i), 'allpts(1,i)', allpts(1,i), 'allpts(2,i)', allpts(2,i)
End Do

!The write_xyy_data subroutine writes the data in the arrays to the specified data file
!Call write_xyy_data ( '/Users/DanielK/Desktop/Fortran/IntRKMethod_'//namein//'.dat', 0, 2, numSteps, allpts, err )

!The write_xyy_plot subroutine can write the script to produce a graph of the data
!Call write_xyy_plots ( '/Users/DanielK/Desktop/Fortran/IntRKMethod.txt', &
!'/Users/DanielK/Desktop/Fortran/IntRKMethod.dat', numSteps, err )

!The run_gnuplot subroutine uses the given script to create a graph of the data
Call run_gnuplot('/Users/DanielK/Desktop/Fortran/IntRKMethod_'//namein//'.txt')

if(err /= 0) then
	Print*, 'Fatal Error: Unable to plot IntRKMethod_'//namein//'.dat'
End if

End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The following subroutines are modified from those provided at the link below:
! http://people.sc.fsu.edu/~jburkardt/f_src/gnufor/gnufor.html
subroutine write_xy_data ( data_filename, n, x, y, ierror )

!*****************************************************************************80
!
!! WRITE_XY_DATA writes X(1:N), Y(1:N) data to a file.
!Note: Editied to write X(0:N), Y(0:N)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the X and Y data
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!

implicit none

integer ( kind = 4 ) n

character ( len = * ) data_filename
integer ( kind = 4 ) file_unit
integer ( kind = 4 ) i
integer ( kind = 4 ) ierror 
integer ( kind = 4 ) ios
real ( kind = 8 ) x(0:n)
real ( kind = 8 ) y(0:n)

ierror = 0

call get_unit ( file_unit )

if ( file_unit == 0 ) then
ierror = 1
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
return
end if

open ( unit = file_unit, file = data_filename, status = 'replace', &
iostat = ios )

if ( ios /= 0 ) then
ierror = 2
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
write ( *, '(a)' ) '  Could not open the output file.'
return
end if

!Edited to start at i=0
do i = 0, n
write ( file_unit, * ) x(i), y(i)
end do

close ( unit = file_unit )

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XY_DATA:'
write ( *, '(a)' ) '  Wrote the GNUPLOT XY data file "' // &
trim ( data_filename ) // '"'

return
end
subroutine write_xy_plot ( command_filename, data_filename, ierror )

!*****************************************************************************80
!
!! WRITE_XY_PLOT writes GNUPLOT commands to plot X(1:N), Y(1:N) data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
implicit none

character ( len = * ) command_filename
character ( len = * ) data_filename
integer ( kind = 4 ) file_unit
integer ( kind = 4 ) i
integer ( kind = 4 ) ierror
integer ( kind = 4 ) ios
!
!  Write the data file.
!
ierror = 0

call get_unit ( file_unit )

if ( file_unit == 0 ) then
ierror = 1
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
return
end if

open ( unit = file_unit, file = command_filename, status = 'replace', &
iostat = ios )

if ( ios /= 0 ) then
ierror = 2
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
write ( *, '(a)' ) '  Could not open the output file.'
return
end if

write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
write ( file_unit, '(a)' ) 'set xlabel "x"'
write ( file_unit, '(a)' ) 'set ylabel "y"'
write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename ) // &
'" using 1:2 with lines'
write ( file_unit, '(a)' ) 'pause -1'
write ( file_unit, '(a)' ) 'q'

close ( unit = file_unit )

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XY_PLOT:'
write ( *, '(a)' ) '  Wrote the GNUPLOT XY plot command file "' // &
trim ( command_filename ) // '"'

return
end

subroutine run_gnuplot ( command_filename )

!*****************************************************************************80
!
!! RUN_GNUPLOT runs GNUPLOT with a given command file.
!
!  Discussion:
!
!    The GNUPLOT program must be available.  To check whether
!    this is so, try typing
!
!      which gnuplot
!
!    If the response is
!
!      gnuplot: command not found
!
!    then you're going to have to make GNUPLOT available.
!
!    You may need to set the environment variable GNUTERM:
!
!      setenv GNUTERM x11
!
!    so that GNUPLOT automatically displays to your X window terminal.
!
!
!    This routine expects that there is a text file containing the appropriate
!    commands to GNUPLOT to display your picture.  There are a number of
!    routines in this package that will do this for simple plotting tasks.
!    Most of them require that you also set up a file of data to be plotted.
!
!    Once this routine invokes GNUPLOT, a graphics window should open
!    up, and the FORTRAN program will pause.  Hitting RETURN should advance
!    to the next picture, or terminate the window at the end, allowing the
!    FORTRAN routine to proceed.
!
!
!    You can look at the data and command files created by the routines.
!    Moreover, you can easily modify the command file to change the options
!    used in GNUPLOT, and then run GNUPLOT interactively, as in:
!
!      gnuplot commands
!
!    In particular, if you want a PostScript version of your graphics files,
!    insert the command "set term postscript" at the beginning of the command
!    file and run gnuplot as follows:
!
!      gnuplot commands > mypicture.ps
!
!    You will also have to hit RETURN once for each plot that is made.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
implicit none

character ( len = 255 ) command
character ( len = * ) command_filename
integer ( kind = 4 ) status
integer ( kind = 4 ) system

call timestamp ( )
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'GNUFOR:'
write ( *, '(a)' ) '  GNUPLOT / FORTRAN90 command interface.'
!
!  Issue a command to the system that will start GNUPLOT, using
!  the file we just wrote as input.
!
!  The "&" will run GNUPLOT in the background, so the FORTRAN program
!  can continue execution independently, and the PERSIST switch tells
!  GNUPLOT that if there are multiple plots, they should each go in
!  a separate window.
!
!  Thanks to Morag Am-Shallem for suggesting these improvements.
!  17 October 2007
!
write ( command, * ) 'gnuplot ' // command_filename

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'GNUFOR:'
write ( *, '(a)' ) '  Issuing the command:"' // trim ( command ) // '".'
write ( *, '(a)' ) ' '
write ( *, '(a)' ) '  Press RETURN to proceed.'

status = system ( trim ( command ) )

if ( status /= 0 ) then
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'GNUFOR - Fatal error!'
write ( *, '(a)' ) '  An error code was returned when the GNUPLOT command'
write ( *, '(a)' ) '  was issued.  Perhaps GNUPLOT is not in your path.'
write ( *, '(a)' ) '  Type "which gnuplot" to check this.'
stop
end if
!
!  Terminate.
!
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'GNUFOR:'
write ( *, '(a)' ) '  Normal end of execution.'

write ( *, '(a)' ) ' '
call timestamp ( )

return
end

subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
implicit none

integer ( kind = 4 ) i
integer ( kind = 4 ) ios
integer ( kind = 4 ) iunit
logical lopen

iunit = 0

do i = 1, 99

if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

inquire ( unit = i, opened = lopen, iostat = ios )

if ( ios == 0 ) then
if ( .not. lopen ) then
iunit = i
return
end if
end if

end if

end do

return
end

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
implicit none

character ( len = 8 ) ampm
integer ( kind = 4 ) d
integer ( kind = 4 ) h
integer ( kind = 4 ) m
integer ( kind = 4 ) mm
character ( len = 9 ), parameter, dimension(12) :: month = (/ &
'January  ', 'February ', 'March    ', 'April    ', &
'May      ', 'June     ', 'July     ', 'August   ', &
'September', 'October  ', 'November ', 'December ' /)
integer ( kind = 4 ) n
integer ( kind = 4 ) s
integer ( kind = 4 ) values(8)
integer ( kind = 4 ) y

call date_and_time ( values = values )

y = values(1)
m = values(2)
d = values(3)
h = values(5)
n = values(6)
s = values(7)
mm = values(8)

if ( h < 12 ) then
ampm = 'AM'
else if ( h == 12 ) then
if ( n == 0 .and. s == 0 ) then
ampm = 'Noon'
else
ampm = 'PM'
end if
else
h = h - 12
if ( h < 12 ) then
ampm = 'PM'
else if ( h == 12 ) then
if ( n == 0 .and. s == 0 ) then
ampm = 'Midnight'
else
ampm = 'AM'
end if
end if
end if

write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

return
end

subroutine write_xyy_data ( data_filename, lda, nrow, ncol, x, ierror )

!*****************************************************************************80
!
!! WRITE_XYY_DATA writes a table of data to a file, for plotting by GNUPLOT.
!
!  Discussion:
!
!    The first column of data is assumed to be the independent variable, X.
!    Separate plots are made of X versus all the other columns of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of X.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the dimensions of X.
!
!    Input, real ( kind = 8 ) X(LDA,NCOL), the NROW by NCOL data to be written.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
implicit none

integer ( kind = 4 ) lda
integer ( kind = 4 ) ncol

character ( len = * ) data_filename
integer ( kind = 4 ) file_unit
integer ( kind = 4 ) i
integer ( kind = 4 ) ierror
integer ( kind = 4 ) ios
integer ( kind = 4 ) nrow
real ( kind = 8 ) x(0:lda,0:ncol)
!Note:Now starts from 0
ierror = 0

call get_unit ( file_unit )

if ( file_unit == 0 ) then
ierror = 1
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
return
end if

open ( unit = file_unit, file = data_filename, status = 'replace', &
iostat = ios )

if ( ios /= 0 ) then
ierror = 2
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
write ( *, '(a)' ) '  Could not open the output file.'
return
end if

!Note: Changed i from starting at 1 to starting at 0
do i = 0, nrow
!Note: Changed from starting at 1 to starting at 0 for x(i,0:ncol)
write ( file_unit, * ) x(i,0:ncol)
end do

close ( unit = file_unit )

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XYY_DATA:'
write ( *, '(a)' ) '  Wrote the GNUPLOT XYY data file "' // &
trim ( data_filename ) // '"'

return
end

subroutine write_xyy_plots ( command_filename, data_filename, &
ncol, ierror )

!*****************************************************************************80
!
!! WRITE_XYY_PLOTS writes GNUPLOT commands to make multiple (X,Y) plots.
!
!  Discussion:
!
!    The first column of data is assumed to be the independent variable, X.
!    Separate plots are made of X versus all the other columns of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns of data.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
implicit none

character ( len = * ) command_filename
character ( len = * ) data_filename
integer ( kind = 4 ) file_unit
integer ( kind = 4 ) i
integer ( kind = 4 ) ierror
integer ( kind = 4 ) ios
integer ( kind = 4 ) ncol
!
!  Write the data file.
!
ierror = 0

call get_unit ( file_unit )

if ( file_unit == 0 ) then
ierror = 1
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
return
end if

open ( unit = file_unit, file = command_filename, status = 'replace', &
iostat = ios )

if ( ios /= 0 ) then
ierror = 2
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
write ( *, '(a)' ) '  Could not open the output file.'
return
end if

write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
write ( file_unit, '(a)' ) 'set xlabel "x"'
write ( file_unit, '(a)' ) 'set ylabel "y"'
do i = 2, ncol
write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename ) // &
'" using ', i, ' with lines'
write ( file_unit, '(a)' ) 'pause -1'
end do
write ( file_unit, '(a)' ) 'q'

close ( unit = file_unit )

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'WRITE_XYY_PLOTS:'
write ( *, '(a)' ) '  Wrote the GNUPLOT XYY plots command file "' // &
trim ( command_filename ) // '"'

return
end
