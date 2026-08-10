# An error-bounded field model of elasmobranch electroreception: canal funnelling, frequency tuning and the effect of body plan in a hammerhead shark and a manta ray

H. Talleb (1)

(1) Sorbonne Université, Université Paris-Saclay, CNRS, CentraleSupélec, Laboratoire de Génie Électrique et Électronique de Paris (GeePs), 91192 Gif-sur-Yvette, France.

Corresponding author: H. Talleb (hakeim.talleb@sorbonne-universite.fr).

Keywords: electroreception, ampullae of Lorenzini, magnetoreception, hammerhead shark, manta ray, cephalofoil, frequency tuning, guaranteed error bounds.

## Abstract

The ampullae of Lorenzini let sharks and rays detect fields of a few nanovolts per centimetre and, through the motional field of swimming, the geomagnetic field. The organ is usually modelled one modality at a time, with a discretisation error estimated rather than bounded. We treat the whole animal as one electroquasistatic conduction problem, in which gel canal, body tissue and seawater are three conductivities, and solve it with a dual discrete geometric method whose complementary energy bounds bracket the answer from both sides. The bracket is a guaranteed interval, so a modelled sensitivity can be compared with a measured behavioural threshold without a numerical caveat. The gel canal raises the certified access conductance of one ampulla from [0.353, 0.367] to [0.576, 0.603] siemens per metre of depth; the intervals are disjoint, so the funnelling gain, guaranteed to lie in [1.57, 1.71], holds for the exact solution. The receptor voltage varies by a factor of six with field direction, so one ampulla is already directional and an array reads field direction. Adding the capacitive sensory epithelium gives a band-pass with corners at 0.2 and 8 hertz peaking near 1.3 hertz, the tuning measured in ampullary afferents. On surface models of a great hammerhead (Sphyrna mokarran) and a manta ray (Mobula birostris), swimming at one metre per second in a fifty microtesla field puts the magnetic modality about a hundred times above threshold. The two body plans then read the same fields differently: the wide cephalofoil converts a uniform navigation field into a large lateral baseline, whereas the compact manta array reads a smaller baseline and a forward-focused prey response. A uniform field and a local prey dipole therefore leave different spatial signatures, separable by geometry rather than amplitude. Read in reverse, the model quantifies the electrical visibility of a body in seawater as a conductivity contrast.

## 1. Introduction

The ampullae of Lorenzini are the most sensitive electroreceptors known [1], [2]. They let elasmobranchs detect the weak bioelectric fields that leak from the prey they hunt and, it is widely held, the Earth's magnetic field through the motional field set up by swimming [3]; the physiology and the ecology of the system have been reviewed recently [4]. Each ampulla is a canal filled with a highly conductive gel, running from a pore in the skin to a sensory epithelium, embedded in resistive body tissue and immersed in seawater, so the organ is, physically, a conduction problem in a piecewise-conductive body.

Three things about that problem have been studied largely apart from one another. The first is the canal. It is a low-resistance path that carries the potential at the pore down to the receptor, so the sensory epithelium reads a difference between a point on the skin and the interior of the animal; its wall resistivity and the electrical properties of the gel have been measured directly [5], [6], [7]. The second is the frequency response. The receptor is band-pass with a peak below ten hertz [8], [9], an origin of which lies in the capacitance of the sensory membrane working against the resistance of the skin and canal wall, and in the oscillatory dynamics of the receptor membrane itself [10]. The third is what the animal does with an array of such receptors, which is where morphology enters. Pore counts and pore distributions vary widely across body plans [11], [12], the wide cephalofoil of hammerhead sharks has long been argued to be an electrosensory adaptation [13], [14], and the peripheral morphology has been proposed as the main determinant of what the electrosensory system can do [15]. That argument is not settled: enhanced binocular overlap [16] and improved manoeuvring [17] are competing explanations for the same head shape. On the magnetic side, behavioural conditioning shows that elasmobranchs detect and can use the geomagnetic field [18], [19], [20], and the ampullary afferents themselves respond to a changing magnetic field [21].

The organ has in fact been modelled repeatedly, and it is worth being precise about what has been done, because it sets what is left to do. Lumped equivalent circuits of the canal and the receptor membrane, including their thermal-noise limits, go back some decades [22], [23]. The field of a prey dipole and its detection by an array have been treated as an electrostatic problem, both to explain strike behaviour [24] and to model prey localisation [25]. The geometry of the canals has been used to derive the directional sensitivity of individual ampullae and the vector coverage of the array [26]. The motional-field route to navigation has been modelled explicitly, including the shorting effect of the surrounding seawater [27]. The most complete treatment to date solves the field over a skate with measured canal geometry, under both a dipole source and a uniform field, and carries it through to the afferent response [28]. This journal has itself published a theory of the compass sense built on the same organ [29].

What none of that supplies is an error bar. The quantities of interest here are receptor voltages of a few hundred nanovolts, set against a measured behavioural threshold of a few nanovolts per centimetre [1]. A conclusion of the form "the modelled signal is above threshold", or "this body plan reads a larger signal than that one", is only as good as the discretisation error of the solve, and a conventional finite-element or finite-difference solution supplies at best an error estimate, not a bound. If the mesh is coarse, as it must be when the canals are millimetric and the animal is metric, the honest statement is that the numerical error is unknown. A convergence study does not close the gap: it shows that a sequence of solutions settles, but never says by how much the last one can still be wrong.

We remove that caveat. We solve the electroquasistatic conduction problem with the dual discrete geometric method [30], and we exploit the fact that a primal solve in the potential and a complementary solve in the current flux bracket the same dissipated energy from above and from below [31]. The interval between them contains the exact value, with no empirical safety factor and no asymptotic assumption. Applied to an ampulla, it brackets the access conductance of the organ, which is the quantity that sets how efficiently the canal funnels an external field to the receptor. The gain is not a tighter number but a different kind of statement: when two bracketed intervals are disjoint, the comparison between them holds for the exact solution and not merely for the computed one.

The contributions are as follows. First, a field model of the ampulla as an electroquasistatic conduction problem with guaranteed two-sided bounds on its access conductance, which is, to our knowledge, the first error-bounded treatment of this organ. Second, the frequency tuning, recovered from the same model by adding the capacitive sensory epithelium. Third, the electric and the magnetic modality quantified on the same engine and referred to the same behavioural threshold, so that their spatial signatures can be compared directly. Fourth, a comparison of two extreme body plans, a hammerhead shark and a manta ray, on photogrammetric surface models, showing that morphology and ampulla placement set the detection geometry. Fifth, the inverse reading of the same model: the electrical visibility of a body in seawater is set by its conductivity contrast, so that a conductivity-matched coating cancels the signature the sense reads while an insulating one increases it.

## 2. Model

### 2.1. The organ as an electroquasistatic conduction problem

At the frequencies of electroreception, below a few tens of hertz, the medium is described by a complex conductivity carrying both conduction and the dielectric effect,

kappa = sigma + j omega epsilon (1)

and the displacement term is negligible against conduction everywhere except across the sensory membrane, where it carries the tuning of section 2.3. The seawater, the body tissue and the ampulla gel are three values of sigma. The gel is taken as conductive as seawater, following the classical description of the canal as a low-resistance core in a resistive body and the direct measurements of the canal wall and of the jelly [5], [6], [7]; the tissue is two orders of magnitude less conductive. Table 1 lists the values, and section 4 discusses what the results owe to this choice.

There are no free charges in the bulk, so the potential obeys current conservation, and the problem is closed by the external stimulus, applied either as a uniform field or as a current dipole in the seawater, and by the receptor and reference electrodes.

### 2.2. Discrete geometric method and guaranteed energy bounds

We discretise with the discrete geometric method, also called the cell method, which writes the field laws directly on a mesh rather than on differential operators. Two interlocked meshes are used: a primal mesh, here a Delaunay tetrahedral complex, and its dual, built on the circumcentres. The scalar potential lives on the primal nodes, its gradient on the primal edges, the current flux on the dual faces, and current conservation on the dual volumes. The differential operators are then incidence matrices, purely topological tables of plus and minus one that carry no metric and no approximation: the gradient is the node-to-edge incidence and the divergence is the dual face-to-volume incidence. All of the geometry and all of the material enter through one object, the constitutive or Hodge operator, which maps an edge quantity to the dual-face quantity through the local conductivity and permittivity, and this is the only place where the discretisation is approximate. Because the topological part is exact, the scheme conserves current cell by cell, which is what makes it amenable to guaranteed error control.

That control is the complementary, or hypercircle, principle. A primal solve takes the potential as the unknown, satisfies conservation exactly and the constitutive law approximately, and over-estimates the dissipated energy. A complementary flux solve takes the current as the unknown, in a divergence-conforming space, satisfies the constitutive law in the dual sense and under-estimates the same energy. The exact value is therefore trapped,

W_lower <= W_exact <= W_upper (2)

with no empirical safety factor. For a voltage-driven problem the access conductance is twice this energy and inherits the same guaranteed bracket. This is the certificate we attach to the funnelling efficiency of the organ, and it is what lets a modelled sensitivity be compared with a behavioural threshold without a discretisation caveat.

### 2.3. The canal, the sensory membrane and the frequency response

The gel canal, almost as conductive as seawater, carries the potential at the pore down to the receptor while the resistive body stays nearly equipotential, so the sensory membrane sees the external potential drop along the canal, read against a deep-body reference. The canal contributes an access resistance. For the lumped frequency model we take it from the analytic resistance of a cylindrical canal, L/(sigma_gel A), which for the dimensions of Table 1 is 15.9 kilohms; this is a different computation from the bracketed conductance of section 2.2, consistent with it in order of magnitude but not identical to it, and the supplementary information sets out the distinction.

The sensory epithelium is a leaky capacitor, a thin interface of surface admittance

Y_m = (sigma_m + j omega epsilon_m)/d_m (3)

which low-passes the response, while the skin and canal wall high-pass it. Their product is the transfer from the external potential to the membrane,

H(omega) = j omega tau_w/(1 + j omega tau_w) . 1/(1 + j omega tau_m) (4)

with tau_w the skin and canal-wall time constant and tau_m the membrane time constant. The response therefore peaks at the geometric mean of the two corners,

f_0 = 1/(2 pi sqrt(tau_w tau_m)) (5)

### 2.4. Electric and magnetic stimuli

A prey is modelled as a bioelectric current dipole in the seawater, whose field falls as the inverse cube of distance. The dipole idealisation and the amplitudes it should carry are supported by direct measurement of the bioelectric fields of marine organisms [32] and by measurement of the charge distribution around a seawater dipole in the behavioural assays themselves [14]. Navigation is modelled by the motional field set up when the animal swims at velocity v through the geomagnetic field B,

E_mot = v x B (6)

which is uniform over the body. The two stimuli are applied with the same engine and read on the same ampulla array, which is what makes their spatial signatures comparable.

### 2.5. Anatomical geometries and meshes (figure 1)

The conductivity field, the body, the seawater and the gel canals are carried by a single unstructured tetrahedral mesh, the primal Delaunay complex, with its circumcentric dual. Two body plans are used, both from photogrammetric surface models: a great hammerhead (Sphyrna mokarran) and a manta ray (Mobula birostris). Each surface model is reoriented to a body frame, scaled to a common body length, and immersed in a seawater box. The ampullae are then placed on the real ventral skin: for the hammerhead they are spread laterally across the ventral cephalofoil, for the manta they are clustered ventrally around the terminal mouth, in both cases by reading the ventral surface height of the actual mesh at each target position and running the canal from that pore into the body. The number of ampullae is reduced to a representative array of nine per animal, so the arrays are anatomically placed but not anatomically complete.

Figure 1 shows the mesh actually used on the hammerhead head: the body immersed in seawater with the ventral array, and a horizontal slab of tetrahedra coloured by material. The canals enter as distinct high-conductivity tubes inside the resistive body rather than as an imposed transfer function, so the funnelling and the directionality come out of the conduction problem itself. Their resolution, however, is not that of section 3.1. The animal is metric and these meshes are built on the whole head, so the canal radius here is comparable to the element size, as the inset of figure 1 shows: the tubes are represented, not resolved. The guaranteed brackets of section 3.1 come from the refined plane model, where the element size is taken below the canal half-width, and the three-dimensional results of sections 3.4 and 3.5 should be read as levels and spatial patterns rather than as certified values. This is the numerical noise referred to in section 4, and a graded mesh on the canals is its remedy.

Table 1. Model parameters.

## 3. Results

### 3.1. Canal funnelling and directional tuning (figure 2)

Figure 2 shows one ampulla under an external field. With the field along the canal, the receptor voltage is 24.2 millivolts per volt per metre of applied field, against 12.3 for the same path without the gel canal: the canal doubles the signal delivered to the receptor.

The certified access conductance makes this quantitative and guaranteed. The problem is plane, so the bracketed quantity is a conductance per unit depth. With the gel canal it is enclosed in [0.576, 0.603] siemens per metre of depth, a half-width of 2.3 per cent, and without it in [0.353, 0.367]. The two intervals are disjoint, so the funnelling gain is itself bracketed, and guaranteed to lie in [1.57, 1.71]. That the canal raises the access conductance is therefore a property of the exact solution, not an artefact of the mesh. The dimensionless gain, being a ratio of two conductances computed on the same geometry, is also the quantity that carries no per-unit-depth caveat, and it is the one swept in the supplementary sensitivity study.

The response is strongly directional. Sweeping the direction of the external field, the receptor voltage varies from 4.9 to 30.8 millivolts per volt per metre, a factor of six between the best and the worst orientation. A single ampulla is therefore already a direction-selective sensor, which is the ingredient the array of section 3.2 uses.

### 3.2. The array as a directional rosette (figure 3)

On a three-dimensional head, an array of ampullae whose canals point along different azimuths gives directional coverage: each ampulla is most sensitive to the field direction aligned with its own canal (figure 3), so the population response encodes field direction. This is the directional rosette of the organ, and it follows from the single-ampulla directionality of section 3.1 without any additional assumption.

### 3.3. Frequency tuning (figure 4)

Adding the capacitive membrane turns the direct-current sensitivity into a band-pass (figure 4). With the values of Table 1, the skin and canal wall give a high-pass corner at 0.20 hertz, the sensory membrane a low-pass corner at 8.0 hertz, and the product peaks at 1.3 hertz, which is within the range of peak sensitivity measured in ampullary afferents [8], [9]. The tuning is therefore not imposed: it follows from the same conduction model once the membrane capacitance is included.

Because the membrane capacitance sets the upper corner (figure 4b), the passband is tunable by that one parameter: a larger capacitance lowers the high edge of the band and narrows it, which gives a concrete physical handle for the interspecific differences in tuning reported in the physiological literature.

### 3.4. Electric and magnetic levels on a hammerhead (figure 5)

Figure 5 puts the two modalities on the hammerhead. Swimming at one metre per second through a fifty microtesla field gives a motional field of fifty microvolts per metre, that is 0.5 microvolts per centimetre, about a hundred times the behavioural threshold of a few nanovolts per centimetre. The magnetic modality is therefore not marginal: on this body plan it is a large signal.

The geometry of the two stimuli differs sharply. Because the receptor voltage tracks lateral position in a uniform field, the wide cephalofoil, which spreads the ventral array over a large lateral span, reads the uniform navigation field as a V-shaped baseline across the array, smallest at the centre and largest toward the wing tips. Across the nine ampullae the navigation baseline runs from 653 nV at the centre of the array to 7.5 and 8.7 microvolts at the two wing tips, a factor of thirteen from centre to tip. A prey dipole, in contrast, produces a localised peak of 214 nV at the ampullae nearest the prey, against an otherwise flat background.

This is the point at which the model says something about the animal rather than about the organ. The uniform navigation field and the localised prey field leave clearly different spatial signatures on the same array, so the two can in principle be separated by the shape of the population response and not only by its amplitude or its frequency content. And the broad span of the hammer is exactly what amplifies the navigation baseline, which is a field-model statement of the electrosensory hypothesis for the cephalofoil.

### 3.5. Hammerhead versus manta (figure 6)

The body plan sets the detection geometry (figure 6), shown here on the two surface models. The hammerhead spreads its ventral array across the wide cephalofoil, over a lateral span four times that of the manta array, and reads the larger navigation baseline: 8.7 microvolts at the wing tip against 5.6 for the manta, and a prey peak of 214 nV against 114. The manta clusters its ampullae ventrally around the terminal mouth, and its prey response is focused at the front.

The size of the effect is worth stating plainly, because it is smaller than the qualitative argument suggests. A fourfold difference in lateral span buys a factor of 1.6 in navigation baseline and 1.9 in prey response, not an order of magnitude. The gain is sublinear because the receptor voltage tracks lateral position in a field that is itself uniform, while the prey field falls as the inverse cube of distance and is therefore set by proximity rather than by span. The body plan shapes the detection geometry, which is the claim; it does not multiply the available signal.

Two extreme body plans therefore read the same fields with different spatial geometries: the spread hammer favours the uniform field, the compact disc favours the local source. Morphology and ampulla placement, and not the receptor alone, shape what the animal can sense, which is the quantitative form of a long-standing argument in the comparative literature [11], [15].

### 3.6. Electrical visibility of a body in seawater (figure 7)

The same model answers an inverse question. A body immersed in seawater perturbs any ambient field because its conductivity differs from the water, and that perturbation is exactly what the ampullae detect. For a body of conductivity sigma_b in seawater sigma_sw the exterior perturbation is that of an induced dipole proportional to the contrast,

alpha = (sigma_b - sigma_sw)/(sigma_b + 2 sigma_sw) (7)

with the volume and the ambient field. Body tissue is less conductive than seawater, so a bare body carries a negative dipole, and this is the field-model counterpart of the classical observation that an object whose conductivity matches the seawater casts no electric shadow and is invisible to the sense. The same conductivity-contrast rule governs the perturbations that weakly electric fish read in active electrolocation [33], and prey are known to exploit it behaviourally, reducing their own electrical conspicuousness when a predator is near [34].

Figure 7 shows the result, computed on the same engine as the dipole-source integral over the mesh. An insulating coating drives the contrast toward its most negative value and increases the signature by a couple of decibels. A conductive coating tuned to the water, on the other hand, makes the coated body a neutral inclusion, for which the exterior dipole vanishes,

sigma_eff(sigma_core, sigma_shell, f) = sigma_sw (8)

where f is the core volume fraction. Sweeping the coating conductivity, the computed signature passes through a deep null, more than thirty decibels below the bare body, at a coating a few times more conductive than seawater. At electroreception frequencies conduction dominates the complex conductivity, omega epsilon being far below sigma in seawater, so the in-band neutralisation is set by the conductivity match alone; the dielectric term refines only the high edge of the band and the transient.

Biologically, this quantifies what makes an object electrically conspicuous to the sense, which is a conductivity contrast and not a size or a metabolic rate. It also gives a counterintuitive practical rule: to lower the electric signature of a body in seawater one needs a conductive, water-matched skin, and an insulating one makes it more visible.

## 4. Discussion

The model unifies what has been studied piecewise. The same engine gives the funnelling of the canal with a guaranteed bracket, the directional selectivity of one ampulla and of the array, the band-pass tuning from the capacitive membrane, the electric and the magnetic sensitivity at known field levels, and the dependence on morphology. The motional-field route to magnetoreception is quantified alongside prey detection, and the two are shown to differ by their spatial signature and not only by their amplitude, which is a testable prediction: an array-level recording should distinguish a uniform stimulus from a local one by the shape of the population response across the cephalofoil.

It is worth stating plainly what is and is not new here. The forward modelling of this organ is not new. Equivalent circuits of the canal and membrane [22], [23], electrostatic treatments of the prey dipole and the array [24], [25], the geometric derivation of ampullary direction vectors [26], and models of the motional-field route to navigation [27] all precede this work, and the field solution of Camperi and colleagues over a skate with measured canal geometry, under both a dipole and a uniform field [28], covers much of the same physical ground. Our contribution sits on top of that literature rather than beside it: those models compute a number, and this one computes an interval that provably contains the number. To that we add the two-body-plan comparison and the inverse, visibility reading of the same solve.

The guaranteed bounds are what make the comparison with behaviour meaningful. A bracket of a couple of per cent on the access conductance is decisive for the questions asked here, because those questions are comparisons: canal against no canal, hammerhead against manta, modelled level against behavioural threshold. In each case the bracket is far narrower than the separation being tested, so the conclusion holds for the exact solution and not only for the computed one. This is a different standard from a convergence study, which shows that a sequence of solutions settles, but never says by how much the last one can still be wrong. It is also, we would argue, the standard this particular biological question calls for, since the whole interest of the organ is that it works within one or two orders of magnitude of the thermal-noise floor [23], where a factor of two in a modelled sensitivity is the difference between a mechanism that works and one that does not.

How far the conclusions travel is set out in the supplementary information, and the short version belongs here. Thirty-five parameter sets were swept: the conductivity of the gel, of the body tissue and of the seawater, the width of the canal, and the mesh size. For each, the guaranteed bracket on the gain was computed and tested against unity. Thirty-two of the thirty-five exclude it, so in each of those cases the canal provably helps. The three that do not are the three that must not: when the gel is given the conductivity of the tissue the canal is no longer a canal, the true gain is exactly one, and the bracket contains one. The method declines to assert an effect where there is none, which is the most direct check available that the bounds are what they claim to be.

Two practical points come out of the same campaign. The bracket does not settle until the element size reaches the canal half-width, and at coarser resolution the brackets at low gel conductivity fail to separate for a purely numerical reason; refining the mesh separates them at every gel conductivity above that of the tissue. And the frequency tuning survives its parameters: perturbing each of the four RC values independently and log-uniformly by up to a factor of 1.5, a band-pass exists in every draw and its peak stays between 1.0 and 2.8 Hz over the central ninety per cent. The tuning conclusion therefore rests on the topology of the circuit rather than on the particular values chosen.

Read in reverse, the model gives a neutralisation rule: the electric signature of a body in seawater is an induced dipole set by its conductivity contrast, an insulating coating increases it, and a conductive water-matched coating cancels it. This is the conduction-limit counterpart of the transformation-optics cloaks of high-frequency electromagnetics, reached here by impedance matching rather than by a metamaterial, and it makes the classical qualitative observation about conductivity-matched objects quantitative.

The limitations are the following, and they bound what should be read into the numbers.

The certified bracket is rigorous for the access conductance, which is a self-energy quantity, and it is that quantity, not the pointwise receptor voltage, that carries a guarantee. The receptor voltages reported in sections 3.4 and 3.5 are computed but not themselves bracketed; certifying them requires a goal-oriented, reciprocity-based formulation, which is the natural next step and is available in principle on the same dual meshes.

The body surfaces are real photogrammetric models, but the electrosensory anatomy on them is not. The ampullae are placed by a parametric rule on the real ventral skin, nine per animal, with canals of a single length and orientation, whereas a real hammerhead carries on the order of a thousand ampullae with a species-specific distribution of canal lengths and directions. The comparison between the two body plans therefore tests the effect of the gross geometry of the array, which is what it is meant to test, and not the effect of the true ampullary map. The absolute per-ampulla values on the coarser three-dimensional meshes also carry numerical noise, for which a graded mesh on the canals is the remedy.

The conductivity of the gel is taken equal to that of seawater. This is the classical idealisation of the canal as a short circuit, and it is the assumption under which the funnelling gain reported in section 3.1 should be read as an upper case: a less conductive gel would lower the access conductance and reduce the gain, without changing its sign, since the canal remains far more conductive than the surrounding tissue. Reconciling the value used here with the published measurements of the jelly and of the canal wall [5], [6], [7] is a refinement the same model accepts unchanged, and one that would be worth making before the absolute conductances are compared with anything measured.

The membrane and skin parameters are representative values chosen to place the corners in the measured range, not values fitted to a species, so the agreement of the peak near 1.3 hertz with the physiological literature should be read as consistency of the mechanism and not as a quantitative fit. The directional tuning of section 3.1 is computed on a stylised two-dimensional cross-section, so the factor of six should be read as the order of the directional selectivity of a single ampulla, not as a species value.

Finally, the visibility result of section 3.6 is the exterior dipole signature of a passive body. A real animal also carries its own bioelectric sources, which a coating does not cancel, and the residual null depth is set by the coating geometry and by the mesh.

## 5. Conclusion

We have presented an error-bounded electroquasistatic model of elasmobranch electroreception which, on a single engine, brackets the funnelling efficiency of the ampulla with guaranteed two-sided bounds, recovers the band-pass frequency tuning from the capacitive sensory membrane, quantifies the electric and the magnetic sensitivity against the behavioural threshold, and shows that the hammerhead and manta body plans read the same fields with different spatial geometries. Because a uniform navigation field and a local prey field produce different spatial signatures on the same array, the model predicts that the two are separable at the population level. Read in reverse, it quantifies the electrical visibility of a body in seawater as a conductivity contrast. Future work will certify the pointwise receptor transfer by reciprocity, refine the ampullary map toward the real distribution of canal lengths and orientations, fit the membrane to measured tuning curves in a named species, and extend the visibility analysis to a body carrying its own bioelectric sources.

## Supplementary information

The supplementary information gives the full discretisation and the construction of the complementary bound, the mesh resolution criterion, the sensitivity and robustness campaign summarised in section 4, and the code repository with a one-command reproduction of every figure.

## Acknowledgements

We thank DigitalLife3D (Jer Bot, digitallife3d.org) for the great hammerhead (Sphyrna mokarran) and manta ray (Mobula birostris) surface models, used under the Creative Commons Attribution-NonCommercial licence (CC BY-NC 4.0).

## CRediT authorship contribution statement

H. Talleb: conceptualisation, methodology, software, formal analysis, investigation, visualisation, writing (original draft), writing (review and editing).

## Funding

This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.

## Declaration of competing interest

The author declares no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Declaration of generative AI and AI-assisted technologies in the writing process

During the preparation of this work the author used Anthropic Claude in order to assist with writing parts of the numerical code and with drafting and language editing of the manuscript. After using this tool, the author reviewed and edited the content as needed and takes full responsibility for the content of the published article. No result, figure, table or numerical value reported here was produced by a generative model: each one is computed by the deposited scripts, which regenerate every figure and every number in this article from the single entry point reproduce_all.py.

## Data availability

All code is publicly available. The repository contains the field engine, the script that produces each figure, the sensitivity campaign of the supplementary information, and a single entry point, reproduce_all.py, that regenerates every figure and every number in this article from scratch. The two surface models used for the body plans are distributed by DigitalLife3D under a Creative Commons Attribution-NonCommercial licence (CC BY-NC 4.0) and are redistributed under the same terms, together with the reorientation and scaling applied to them. Repository and archived release are given in section S6 of the supplementary information.

## References

[1] A. J. Kalmijn, "The electric sense of sharks and rays," J. Exp. Biol., vol. 55, no. 2, pp. 371-383, 1971. doi:10.1242/jeb.55.2.371

[2] R. W. Murray, "The response of the ampullae of Lorenzini of elasmobranchs to electrical stimulation," J. Exp. Biol., vol. 39, no. 1, pp. 119-128, 1962. doi:10.1242/jeb.39.1.119

[3] Ad. J. Kalmijn, "Electric and magnetic field detection in elasmobranch fishes," Science, vol. 218, no. 4575, pp. 916-918, 1982. doi:10.1126/science.7134985

[4] K. C. Newton, A. B. Gill, and S. M. Kajiura, "Electroreception in marine fishes: chondrichthyans," J. Fish Biol., vol. 95, no. 1, pp. 135-154, 2019. doi:10.1111/jfb.14068

[5] B. Waltman, "Electrical properties and fine structure of the ampullary canals of Lorenzini," Acta Physiol. Scand. Suppl., vol. 264, pp. 1-60, 1966.

[6] B. R. Brown, J. C. Hutchison, M. E. Hughes, D. R. Kellogg, and R. W. Murray, "Electrical characterization of gel collected from shark electrosensors," Phys. Rev. E, vol. 65, no. 6, art. 061903, 2002. doi:10.1103/PhysRevE.65.061903

[7] E. E. Josberger, P. Hassanzadeh, Y. Deng, J. Sohn, M. J. Rego, C. T. Amemiya, and M. Rolandi, "Proton conductivity in ampullae of Lorenzini jelly," Sci. Adv., vol. 2, no. 5, art. e1600112, 2016. doi:10.1126/sciadv.1600112

[8] J. C. Montgomery, "Frequency response characteristics of primary and secondary neurons in the electrosensory system of the thornback ray," Comp. Biochem. Physiol. A, vol. 79, no. 1, pp. 189-195, 1984. doi:10.1016/0300-9629(84)90731-X

[9] T. C. Tricas and J. G. New, "Sensitivity and response dynamics of elasmobranch electrosensory primary afferent neurons to near threshold fields," J. Comp. Physiol. A, vol. 182, no. 1, pp. 89-101, 1998. doi:10.1007/s003590050161

[10] W. T. Clusin and M. V. L. Bennett, "The oscillatory responses of skate electroreceptors to small voltage stimuli," J. Gen. Physiol., vol. 73, no. 6, pp. 685-702, 1979. doi:10.1085/jgp.73.6.685

[11] S. M. Kajiura, "Head morphology and electrosensory pore distribution of carcharhinid and sphyrnid sharks," Environ. Biol. Fishes, vol. 61, no. 2, pp. 125-133, 2001. doi:10.1023/A:1011028312787

[12] R. M. Kempster, I. D. McCarthy, and S. P. Collin, "Phylogenetic and ecological factors influencing the number and distribution of electroreceptors in elasmobranchs," J. Fish Biol., vol. 80, no. 5, pp. 2055-2088, 2012. doi:10.1111/j.1095-8649.2011.03214.x

[13] S. M. Kajiura and K. N. Holland, "Electroreception in juvenile scalloped hammerhead and sandbar sharks," J. Exp. Biol., vol. 205, no. 23, pp. 3609-3621, 2002. doi:10.1242/jeb.205.23.3609

[14] S. M. Kajiura and T. P. Fitzgerald, "Response of juvenile scalloped hammerhead sharks to electric stimuli," Zoology, vol. 112, no. 4, pp. 241-250, 2009. doi:10.1016/j.zool.2008.07.001

[15] T. C. Tricas, "The neuroecology of the elasmobranch electrosensory world: why peripheral morphology shapes behavior," Environ. Biol. Fishes, vol. 60, no. 1-3, pp. 77-92, 2001. doi:10.1023/A:1007684404669

[16] D. M. McComb, T. C. Tricas, and S. M. Kajiura, "Enhanced visual fields in hammerhead sharks," J. Exp. Biol., vol. 212, no. 24, pp. 4010-4018, 2009. doi:10.1242/jeb.032615

[17] S. M. Kajiura, J. B. Forni, and A. P. Summers, "Maneuvering in juvenile carcharhinid and sphyrnid sharks: the role of the hammerhead shark cephalofoil," Zoology, vol. 106, no. 1, pp. 19-28, 2003. doi:10.1078/0944-2006-00086

[18] C. G. Meyer, K. N. Holland, and Y. P. Papastamatiou, "Sharks can detect changes in the geomagnetic field," J. R. Soc. Interface, vol. 2, no. 2, pp. 129-130, 2005. doi:10.1098/rsif.2004.0021

[19] K. C. Newton and S. M. Kajiura, "Magnetic field discrimination, learning, and memory in the yellow stingray (Urobatis jamaicensis)," Anim. Cogn., vol. 20, no. 4, pp. 603-614, 2017. doi:10.1007/s10071-017-1084-8

[20] B. A. Keller, N. F. Putman, R. D. Grubbs, D. S. Portnoy, and T. P. Murphy, "Map-like use of Earth's magnetic field in sharks," Curr. Biol., vol. 31, no. 13, pp. 2881-2886.e3, 2021. doi:10.1016/j.cub.2021.03.103

[21] G. N. Akoev, O. B. Ilyinsky, and P. M. Zadan, "Responses of electroreceptors (ampullae of Lorenzini) of skates to electric and magnetic fields," J. Comp. Physiol. A, vol. 106, no. 2, pp. 127-136, 1976. doi:10.1007/BF00620494

[22] H. M. Fishman, "Appendix: the biophysics of electroreception in ampullary organs of elasmobranch fishes," in Cell Physiology Source Book, 2nd ed. San Diego, CA, USA: Academic Press, 2001, pp. 857-861. doi:10.1016/B978-0-08-057455-4.50069-4

[23] R. K. Adair, R. D. Astumian, and J. C. Weaver, "Detection of weak electric fields by sharks, rays, and skates," Chaos, vol. 8, no. 3, pp. 576-587, 1998. doi:10.1063/1.166339

[24] B. R. Brown, "Modeling an electrosensory landscape: behavioral and morphological optimization in elasmobranch prey capture," J. Exp. Biol., vol. 205, no. 7, pp. 999-1007, 2002. doi:10.1242/jeb.205.7.999

[25] D. Kim, "Prey detection mechanism of elasmobranchs," BioSystems, vol. 87, no. 2-3, pp. 322-331, 2007. doi:10.1016/j.biosystems.2006.09.029

[26] A. C. Rivera-Vicente, J. Sewell, and T. C. Tricas, "Electrosensitive spatial vectors in elasmobranch fishes: implications for source localization," PLoS ONE, vol. 6, no. 1, art. e16008, 2011. doi:10.1371/journal.pone.0016008

[27] T. C. A. Molteno and W. L. Kennedy, "Navigation by induction-based magnetoreception in elasmobranch fishes," J. Biophys., vol. 2009, art. 380976, pp. 1-6, 2009. doi:10.1155/2009/380976

[28] M. Camperi, T. C. Tricas, and B. R. Brown, "From morphology to neural information: the electric sense of the skate," PLoS Comput. Biol., vol. 3, no. 6, art. e113, 2007. doi:10.1371/journal.pcbi.0030113

[29] M. G. Paulin, "Electroreception and the compass sense of sharks," J. Theor. Biol., vol. 174, no. 3, pp. 325-339, 1995. doi:10.1006/jtbi.1995.0102

[30] E. Tonti, "A direct discrete formulation of field laws: the cell method," CMES Comput. Model. Eng. Sci., vol. 2, no. 2, pp. 237-258, 2001. doi:10.3970/cmes.2001.002.237

[31] J. L. Synge, The Hypercircle in Mathematical Physics: A Method for the Approximate Solution of Boundary Value Problems. Cambridge, U.K.: Cambridge Univ. Press, 1957.

[32] C. N. Bedore and S. M. Kajiura, "Bioelectric fields of marine organisms: voltage and frequency contributions to detectability by electroreceptive predators," Physiol. Biochem. Zool., vol. 86, no. 3, pp. 298-311, 2013. doi:10.1086/669973


[33] B. Rasnow, "The effects of simple objects on the electric field of Apteronotus," J. Comp. Physiol. A, vol. 178, no. 3, pp. 397-411, 1996. doi:10.1007/BF00193977

[34] C. N. Bedore, S. M. Kajiura, and S. Johnsen, "Freezing behaviour facilitates bioelectric crypsis in cuttlefish faced with predation risk," Proc. R. Soc. B, vol. 282, no. 1820, art. 20151886, 2015. doi:10.1098/rspb.2015.1886
