      subroutine set_loop_limits

      use mainvar_module

	if(dim.eq.2)then
	    js1=1
	    je1=endy

	    js2=2
	    je2=endy-1

	    js3=3
	    je3=endy-2

	    js5=5
		je5=endy-4

	    js6=6
		je6=endy-5

		js7=7
		je7=endy-6
	else
		js1=1+overlap
		je1=endy-overlap

		js2=1+overlap
		je2=endy-overlap

		js3=1+overlap
		je3=endy-overlap

		js5=1+overlap
		je5=endy-overlap

		js6=1+overlap
		je6=endy-overlap

		js7=1+overlap
		je7=endy-overlap
	endif




      return

      end

