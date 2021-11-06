// Nbody simulation to create cluster of sound sources including position and doppler effect.
// The sound source is given through the synth
//
// v0.4 2016.08.28 : sleon

(

~nsources = 30 ;
~nInternalIteration = 300 ; // internal loop before changing the sound parameters
~internalDt = 5e-5        ; // internal dt for the nbody integration
~soft2 = 0.05            ; // softening ...

~synthId   = 5            ; // Synth to be used
~freqBasis = 100.0          ; // frequency fundamental for some of the synth

~hscale    = 1.0             ; // lenght scale for the gravitational wall ...
~alpha     = -1.0          ; // power law index for the wall .. (not used with exponentional)
~dwall     = 2.0           ; // strenght of Wall


~csound            = 30.0           ;     // sound velocity
~massOneBody       = 2.0       ;       // mass of one body
~posStart          = [0.0,10, 0.0] ;    // initial position
~positionStartDisp = 4.0    ;            // initial position dispersion
~velStart          = [0.0, -2.0 ,0.0] ;     // initial velocity
~velocityStartDisp = 1.5   ;            // initial velocity dispersion


~posx = Array.new(~nsources) ;
~posy = Array.new(~nsources) ;
~posz = Array.new(~nsources) ;

~velx = Array.new(~nsources) ;
~vely = Array.new(~nsources) ;
~velz = Array.new(~nsources) ;

~forcex = Array.new(~nsources) ;
~forcey = Array.new(~nsources) ;
~forcez = Array.new(~nsources) ;

~mass = Array.new(~nsources) ;
~neighborIndex = Array.new(~nsources) ;

~nsynth = Array.new(~nsources) ;

///Synth ..
///


SynthDef(\nbodysound1, {arg inBuf = 8 , pos = #[0.0,10.0,0.0], vel = #[0.0,0.0,0.0] ,force = #[0.0, 0.0, 0.0] ;

    var  env ;
	var amp = 1.0 ;
	var xout, yout ;
	var vproj, r2  ;
	var direction = 1.0 ;
	var dopplerShift = 1.0 ;

	var knorm = 2.0 ;


	r2 = (pos[0]*pos[0]) + (pos[1]*pos[1]) + (pos[2]*pos[2]) ;
	d = r2.sqrt() ;
	amp = knorm / (d*d) ;
	vproj = ((pos[0]*vel[0])+(pos[1]*vel[1])+(pos[2]*vel[2])) / d  ;
	dopplerShift =   ~csound / (~csound + vproj)    ;
	direction    =   atan2(pos[0],pos[1]) / pi ;


	#x, y = FreqShift.ar(In.ar(inBuf,2), 100.0 * dopplerShift) ;

	#xout, yout = Rotate2.ar(x, y, direction) ;

	Out.ar(0,[xout,yout]) ;
}).add ;

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

SynthDef(\nbodysound2, {arg freqFundamental = 55, freqMul = 1, inBuf = 8 , pos = #[0.0,10.0,0.0], vel = #[0.0,0.0,0.0], force = #[0.0, 0.0, 0.0] ;

    var  env, freq0 ;
	var amp = 1.0 ;
	var xout, yout ;
	var vproj, r2  ;
	var direction = 1.0 ;
	var dopplerShift = 1.0 ;
	var knorm = 2.0 ;

	freq0 = freqFundamental * freqMul ;

	r2 = (pos[0]*pos[0]) + (pos[1]*pos[1]) + (pos[2]*pos[2]) ;
	d = r2.sqrt() ;
	amp = knorm / (d*d) ;
	vproj = ((pos[0]*vel[0])+(pos[1]*vel[1])+(pos[2]*vel[2])) / d  ;
	dopplerShift =   ~csound / (~csound + vproj)    ;
	direction    =   atan2(pos[0],pos[1]) / pi ;

	f = SinOsc.ar(freq0*dopplerShift, amp) ;

	#x, y = [f,f] ;
	#xout, yout = Rotate2.ar(x, y, direction) ;

	Out.ar(0,[xout,yout]) ;
}).add ;


///////////////////////////////////////////////////
///////////////////////////////////////////////////

SynthDef(\nbodysound3, {arg freqFundamental = 55, freqMul = 1, inBuf = 8 , pos = #[0.0,10.0,0.0], vel = #[0.0,0.0,0.0], force = #[0.0, 0.0, 0.0] ;

    var  env, freq0 ;
	var amp = 1.0 ;
	var xout, yout ;
	var vproj, r2  ;
	var direction = 1.0 ;
	var dopplerShift = 1.0 ;
	var knorm = 2.0 ;

	freq0 = freqFundamental * freqMul ;

	r2 = (pos[0]*pos[0]) + (pos[1]*pos[1]) + (pos[2]*pos[2]) ;
	d = r2.sqrt() ;
	amp = knorm / (d*d) ;
	vproj = ((pos[0]*vel[0])+(pos[1]*vel[1])+(pos[2]*vel[2])) / d  ;
	dopplerShift =   ~csound / (~csound + vproj)    ;
	direction    =   atan2(pos[0],pos[1]) / pi ;

	f = LFSaw.ar(freq0*dopplerShift, amp) ;

	#x, y = [f,f] ;
	#xout, yout = Rotate2.ar(x, y, direction) ;

	Out.ar(0,[xout,yout]) ;
}).add ;

///////////////////////////////////////////////////////////////////////////////7
///////////////////////////////////////////////////////////////////////////////
//////////// F Modulation
SynthDef(\nbodysound4, {arg freqFundamental = 55, freqMul = 1, inBuf = 8 , pos = #[0.0,10.0,0.0], vel = #[0.0,0.0,0.0], force = #[0.0, 0.0, 0.0] ;

    var  env, freq0 ;
	var amp = 1.0 ;
	var xout, yout ;
	var vproj, r2  ;
	var mod ;
	var carPartial, modPartial ;
	var direction = 1.0 ;
	var dopplerShift = 1.0 ;
	var knorm = 2.0 ;


	r2 = (pos[0]*pos[0]) + (pos[1]*pos[1]) + (pos[2]*pos[2]) ;
	d = r2.sqrt() ;
	amp = knorm / (d*d) ;
	vproj = ((pos[0]*vel[0])+(pos[1]*vel[1])+(pos[2]*vel[2])) / d  ;
	dopplerShift =   ~csound / (~csound + vproj)    ;
	direction    =   atan2(pos[0],pos[1]) / pi ;

	freq0 = freqFundamental * freqMul ;
	carPartial = Rand(0.1,5.0) ;
    modPartial = Rand(0.1,5.0) ;

    mod = SinOsc.ar(
        freq0 * modPartial * dopplerShift,
        0,
        freq0 * freqMul * LFNoise1.kr(5.reciprocal).abs
    );

	f = SinOsc.ar(
		(freq0 * carPartial * dopplerShift) + mod + LFNoise2.kr(0.5,10, -20) , 0, amp );

	#x, y = [f,f] ;
	#xout, yout = Rotate2.ar(x, y, direction) ;

	Out.ar(0,[xout,yout]) ;
}).add ;


///////////////////////////////////////
///////////////////////////////////////

SynthDef(\nbodysound5, {arg freqFundamental = 55, freqMul = 1, inBuf = 8 , pos = #[0.0,10.0,0.0], vel = #[0.0,0.0,0.0], force = #[0.0, 0.0, 0.0] ;

    var  env, freq0 ;
	var amp = 1.0 ;
	var xout, yout ;
	var vproj, r2 , ft ;
	var direction = 1.0 ;
	var dopplerShift = 1.0 ;

	var knorm = 1.0 ;
	var carPar  ;

	freq0 = freqFundamental * freqMul ;

	r2 = (pos[0]*pos[0]) + (pos[1]*pos[1]) + (pos[2]*pos[2]) + ~soft2 ;

	d = r2.sqrt() ;
	amp = knorm / (d*d) ;
	vproj = ((pos[0]*vel[0])+(pos[1]*vel[1])+(pos[2]*vel[2])) / d  ;
	dopplerShift =   ~csound / (~csound + vproj)    ;
	direction    =   atan2(pos[0],pos[1]) / pi ;


	ft = (force[0]*force[0]) + (force[1]*force[1]) + (force[2]*force[2]) ;
		carPar = max(0.1, ft) ;


	s = SinOsc.ar(freq0 * dopplerShift, 0, freq0 * dopplerShift * LFNoise1.kr(5.reciprocal).abs) ;
	u = SinOsc.ar( (freq0 * 0.6 * dopplerShift) + s, 0, amp /3) ;
	v = SinOsc.ar( (freq0 * 0.3 * dopplerShift) + s,0,amp/4) ;
	w = SinOsc.ar( (freq0 * 0.9 * dopplerShift) + s,0,amp/5) ;
	x = SinOsc.ar( (freq0 * 1.2 * dopplerShift) + s,0,amp/6) ;



	f = u + v + w + x + CombL.ar( u+v+w+x, 0.2, ft / 100.0, ft /10.0 ) + CombL.ar( u+v+w+x, 0.2, ft / 80.0, ft / 8.0 );



	#x, y = [f,f] ;
	#xout, yout = Rotate2.ar(x, y, direction) ;

	Out.ar(0,[xout,yout]) ;
}).add ;

////////////////////////////////////
/// Initial conditions
/// Spherical cluster

~initialConditions = {arg	massBody, initPos, posDisp, initVel, velDisp ;

	"## initial conditions...".postln ;

	~nsources.postln ;
	a = Pgauss(initPos[0], posDisp,  inf);
	~posx = a.asStream.nextN(~nsources);
	a = Pgauss(initPos[1], posDisp,  inf);
	~posy = a.asStream.nextN(~nsources);
	a = Pgauss(initPos[2], posDisp,  inf);
	~posz = a.asStream.nextN(~nsources);

	a = Pgauss(initVel[0], velDisp,  inf);
	~velx = a.asStream.nextN(~nsources);
	a = Pgauss(initVel[1], velDisp,  inf);
	~vely = a.asStream.nextN(~nsources);
	a = Pgauss(initVel[2], velDisp,  inf);
	~velz = a.asStream.nextN(~nsources);


	//  .. initialiasing ..
	~mass = a.asStream.nextN(~nsources);
	~forcex = a.asStream.nextN(~nsources);
	~forcey = a.asStream.nextN(~nsources);
	~forcez = a.asStream.nextN(~nsources);
	~nsynth = a.asStream.nextN(~nsources);
	~neighborIndex =  a.asStream.nextN(~nsources);

	for (0, ~nsources-1, { arg i;
		~mass[i] = massBody ;
	}) ;

	"## Ending initial conditions ..".postln ;
} ;


/// Force computation for all bodies
///

~forceIteration = {

	var fx, fy, fz ;
	var r, r2, r3 ;
	var xx, yy, zz , yi ;

//	"### Force iteration ...".postln ;
	for (0, ~nsources-1, {arg i ; ~neighborIndex[i] = i}) ;

	for (0, ~nsources-1, { arg i;
		fx = 0.0 ;
		fy = 0.0 ;
		fz = 0.0 ;

		n = ~neighborIndex.rotate(-1*i-1) ;

		for (0,~nsources-2, {arg j;

			xx = ~posx[i]-~posx[n[j]] ;
			yy = ~posy[i]-~posy[n[j]] ;
			zz = ~posz[i]-~posz[n[j]] ;

			r2 = (xx*xx) + (yy*yy) + (zz*zz) + ~soft2;

			r = r2.sqrt() ;
			r3 = r * r *r ;

			fx = fx - ((~posx[i]-~posx[n[j]]) * (~mass[i]* ~mass[n[j]]) / r3) ;
			fy = fy - ((~posy[i]-~posy[n[j]]) * (~mass[i]* ~mass[n[j]]) / r3) ;
			fz = fz - ((~posz[i]-~posz[n[j]]) * (~mass[i]* ~mass[n[j]]) / r3) ;



		}) ;

		// Wall computation

		yi = ~posy[i] ;

		// fy = fy - ((sign(yi) * (~dwall *  exp( ~alpha * log(((abs(yi)+ ~soft2) / ~hscale)))))) ;

		fy = fy - ((sign(yi) * (~dwall *  exp(-1 * (abs(yi)/ ~hscale))))) ;


		~forcex[i] = fx ;
		~forcey[i] = fy ;
		~forcez[i] = fz ;
	});
} ;

/// Nbody integrator and sound emission
///
///

~main = Routine({

			var timeToWait = 0.5 ;

	        // Intitial Conditions
	        //
			~initialConditions.value(massBody:~massOneBody , initPos:~posStart, posDisp: ~positionStartDisp, initVel:~velStart,velDisp:~velocityStartDisp) ;

			// initiate synth...
	        "###".postln ;
	        "### Creation of synths ...".postln ;
			~nsources.do({arg i ;
		         switch(~synthId,
			1, {~nsynth[i] = Synth.new("nbodysound1")},
			2, {~nsynth[i] = Synth.new("nbodysound2", [\freqFundamental, ~freqBasis, \freqMul, i])},
			3, {~nsynth[i] = Synth.new("nbodysound3", [\freqFundamental, ~freqBasis, \freqMul, i])},
			4, {~nsynth[i] = Synth.new("nbodysound4", [\freqFundamental, ~freqBasis, \freqMul, i])},
			5, {~nsynth[i] = Synth.new("nbodysound5", [\freqFundamental, ~freqBasis, \freqMul, i])}
		       )
			}) ;

			inf.do({arg j ;

				"## Internal computation of forces ...".postln ;

				// internal body iteration
				~nInternalIteration.do({

			      ~forceIteration.value ;

					~nsources.do({arg i;

				        ~velx[i] = ~velx[i] + ((~internalDt * ~forcex[i]) / ~mass[i]) ;
				        ~vely[i] = ~vely[i] + ((~internalDt * ~forcey[i]) / ~mass[i]) ;
				        ~velz[i] = ~velz[i] + ((~internalDt * ~forcez[i]) / ~mass[i]) ;

				        ~posx[i] = ~posx[i] + (~internalDt * ~velx[i]) ;
				        ~posy[i] = ~posy[i] + (~internalDt * ~vely[i]) ;
				        ~posz[i] = ~posz[i] + (~internalDt * ~velz[i]) ;

					});
		         });

				// update the synths ...

		        "### Update the synths ..".postln ;

		        ~posx.postln;
		        ~posy.postln;
		        ~posz.postln;
		        "---".postln ;
		        ~velz.postln;
				~vely.postln;
				~velz.postln;

				~nsources.do({arg i;
			~nsynth[i].set("pos", [~posx[i],~posy[i],~posz[i]], "vel",[~velx[i],~vely[i],~velz[i]] , "force", [~forcex[i],~forcey[i],~forcez[i]]) ;
					}) ;

		        // timeToWait = 1.0.rand() ;
				timeToWait.wait;

			});

}) ;


~main.play ;

)
