# AMLACAL

This is the wiki page of number counting project affiliated of ALMACAL



### Introduction

1. The evolution of the number counts, whether we have resolved the all the extragalactic background
2. The nature of the faint SMG, counterparts of normal star-forming galaxies at high-z
3. the clutering of the faint SMG



### Advantages:

1. **Well-caibrated data**. Because the calibrator, we know whether the field is calibrated of not. Which means we can throw out the data corrupted by some unknown reasons
2. **Small cosmic covariance**. The calibrator are distributed in all the sky, which means the comic variance can be small comparing to the several deep field
3. **Self-calibration** to improve image fidelity 
4. **High-resolution**, reducing the uncertainty introduced by the unresolved multiple sources



### Disadvantage:

1. The calibrator are most blazers, which means that their flux may change from time to time. 
   - **One way to do is to subtract the point source in uv space in the each observation**
   - Directly imaging to get the average 
   - Try to model the variability of the blazer and then scale their flux?

2. The bias introduced by the Calibrator
   - Selection bias around blazer
   - Possible contamination from the resolve structure of the blazer 



[workflow](https://mermaid-js.github.io/mermaid-live-editor/#/edit/eyJjb2RlIjoiZmxvd2NoYXJ0IFREXG4gICAgQUxNQUNBTFtBTE1BQ0FMIG51bWJlciBjb3VudHNdIC0tPiBEYXRhWyhEYXRhKV1cbiAgICBBTE1BQ0FMW0FMTUFDQUwgbnVtYmVyIGNvdW50c10gLS0-IFNpbXVsYXRpb25bW1NpbXVsYXRpb25dXVxuIFxuICAgIHN1YmdyYXBoIFN0ZXAxIFtEYXRhIHByZS1wcm9jZXNzaW5nXVxuICAgIHV2ZGF0YSAtLT4gczFfMShjYWxpYnJhdGUgdGhlIGNhaWxpYnJhdG9yKVxuICAgIHV2ZGF0YSAtLT4gczFfMihzcGxpdCB0aGUgY2FsaWJyYXRvcilcbiAgICB1dmRhdGEgLS0-IHMxXzMocG9pbnQgc291cmNlIHN1YnRyYWN0aW9uKVxuICAgIHNwdyAtLT58Pz8_fCBhbGlnbm1lbnQoYWxpZ25lbWVudClcbiAgICBcbiAgICBlbmRcbiAgICBcbiAgICBzdWJncmFwaCBTdGVwMiBbSW1hZ2luZ11cbiAgICBpbWFnaW5nIC0tPiBzMl8xW2luZGl2aWR1YWwgaW1hZ2VzXVxuICAgIGltYWdpbmcgLS0-IHMyXzJbY29uYmluYXRpb25dXG4gICAgIHNvdXJjZV9leHRyYWN0aW9uW1NvdXJjZSBleHRyYWN0aW9uXSAtLT4gcmVzb2x1dGlvbihbUmVzb2x1dGlvbl0pXG4gICAgIHNvdXJjZV9leHRyYWN0aW9uW1NvdXJjZSBleHRyYWN0aW9uXSAtLT4gZWZmZWN0aXZlX2FyZWEoW0VmZmVjdGl2ZSBBcmVhXSlcbiAgICAgc291cmNlX2V4dHJhY3Rpb25bU291cmNlIGV4dHJhY3Rpb25dIC0tPiBzZXh0cmFjdG9yKFtzZXh0cmFjdG9yXSlcbiAgICBlbmRcblxuICAgIHN1YmdyYXBoIFN0ZXAzIFtOdW1iZXIgQ291bnRpbmddXG4gICAgc3VyZW5lc3NbU3VyZW5lc3NdIC0tPiBtdWx0aWJhbmQoW011bHRpYmFuZCBjb3VudGVycGFydF0pXG4gICAgZW5kXG4gIFxuXG4gICAgRGF0YSA9PT4gU3RlcDFcbiAgICBTdGVwMSA9PT4gU3RlcDJcbiAgICBTdGVwMiA9PT4gU3RlcDNcblxuICAgIFNpbXVsYXRpb24gLS0-IFRlc3RzXG4gICAgc3ViZ3JhcGggVGVzdHMgW1Rlc3RpbmddXG4gICAgczFbL0NvbXBsZXRlbmVzc1xcXSAtLS0gczJbL0ZsdXggYm9vc3RpbmdcXF0gLS0tfD8_P3wgczNbL1NlbGVjdGlvbiBiaWFzXFxdXG4gICAgZW5kXG4gICAgVGVzdHMgLS4tPiBTdGVwM1xuICAgIFN0ZXAzIC0tPiBjb21wYXJpc29tW0NvbXBhcmluZyB3aXRoIG1vZGVsc10iLCJtZXJtYWlkIjp7fSwidXBkYXRlRWRpdG9yIjpmYWxzZX0)

![hello](https://mermaid.ink/svg/eyJjb2RlIjoiZmxvd2NoYXJ0IFREXG4gICAgQUxNQUNBTFtBTE1BQ0FMIG51bWJlciBjb3VudHNdIC0tPiBEYXRhWyhEYXRhKV1cbiAgICBBTE1BQ0FMW0FMTUFDQUwgbnVtYmVyIGNvdW50c10gLS0-IFNpbXVsYXRpb25bW1NpbXVsYXRpb25dXVxuIFxuICAgIHN1YmdyYXBoIFN0ZXAxIFtEYXRhIHByZS1wcm9jZXNzaW5nXVxuICAgIHV2ZGF0YSAtLT4gczFfMShjYWxpYnJhdGUgdGhlIGNhaWxpYnJhdG9yKVxuICAgIHV2ZGF0YSAtLT4gczFfMihzcGxpdCB0aGUgY2FsaWJyYXRvcilcbiAgICB1dmRhdGEgLS0-IHMxXzMocG9pbnQgc291cmNlIHN1YnRyYWN0aW9uKVxuICAgIHNwdyAtLT58Pz8_fCBhbGlnbm1lbnQoYWxpZ25lbWVudClcbiAgICBcbiAgICBlbmRcbiAgICBcbiAgICBzdWJncmFwaCBTdGVwMiBbSW1hZ2luZ11cbiAgICBpbWFnaW5nIC0tPiBzMl8xW2luZGl2aWR1YWwgaW1hZ2VzXVxuICAgIGltYWdpbmcgLS0-IHMyXzJbY29uYmluYXRpb25dXG4gICAgIHNvdXJjZV9leHRyYWN0aW9uW1NvdXJjZSBleHRyYWN0aW9uXSAtLT4gcmVzb2x1dGlvbihbUmVzb2x1dGlvbl0pXG4gICAgIHNvdXJjZV9leHRyYWN0aW9uW1NvdXJjZSBleHRyYWN0aW9uXSAtLT4gZWZmZWN0aXZlX2FyZWEoW0VmZmVjdGl2ZSBBcmVhXSlcbiAgICAgc291cmNlX2V4dHJhY3Rpb25bU291cmNlIGV4dHJhY3Rpb25dIC0tPiBzZXh0cmFjdG9yKFtzZXh0cmFjdG9yXSlcbiAgICBlbmRcblxuICAgIHN1YmdyYXBoIFN0ZXAzIFtOdW1iZXIgQ291bnRpbmddXG4gICAgc3VyZW5lc3NbU3VyZW5lc3NdIC0tPiBtdWx0aWJhbmQoW011bHRpYmFuZCBjb3VudGVycGFydF0pXG4gICAgZW5kXG4gIFxuXG4gICAgRGF0YSA9PT4gU3RlcDFcbiAgICBTdGVwMSA9PT4gU3RlcDJcbiAgICBTdGVwMiA9PT4gU3RlcDNcblxuICAgIFNpbXVsYXRpb24gLS0-IFRlc3RzXG4gICAgc3ViZ3JhcGggVGVzdHMgW1Rlc3RpbmddXG4gICAgczFbL0NvbXBsZXRlbmVzc1xcXSAtLS0gczJbL0ZsdXggYm9vc3RpbmdcXF0gLS0tfD8_P3wgczNbL1NlbGVjdGlvbiBiaWFzXFxdXG4gICAgZW5kXG4gICAgVGVzdHMgLS4tPiBTdGVwM1xuICAgIFN0ZXAzIC0tPiBjb21wYXJpc29tW0NvbXBhcmluZyB3aXRoIG1vZGVsc10iLCJtZXJtYWlkIjp7InRoZW1lIjoiZGVmYXVsdCJ9LCJ1cGRhdGVFZGl0b3IiOmZhbHNlfQ)



## Recordings of every day:

**20120-12-07**

The calibrator with Jet: J0635-7156

1. Careful manual calibration, try to find the difference

   - 2016.1.00203.S   B6   member.uid___A001_X87c_X8bf >> 

   - filename on HLTau: uid___A002_Xbf4033_X9d8.ms.split.cal.J0635-7516_B6

2. 10good vs 10bad, calibrator: J1717-3342

----------------------

**2020-12-08**

possibe difference between good residual and bad residuals (based on J0217-0820):

- [x] ~~integration time~~
- [x] ~~bandpass or gain calibrator~~
- [x] **uv coverage**, the good ones tend to have higher resolution
- [x] ~~spw window coverage~~. (maybe update the spw_stat to handle multiple vis)
- [x] ~~elevation angle~~

**thoughts**: 

1. ~~the non-detection could also be used to constrain the number density of the SMG, assuming a model~~, it is already considered when we calculate the effective area

-----------------------

**2020-12-09**

Plans for this week and maybe the following weeks:

- [ ] Better analysis the field of J0217-0820
  - [ ] test tclean with multi-threshold method
  - [ ] auto source detection
  - [ ] multi-band auto source confirmation
  - [ ] made the completeness and flux bootsting simulation
- [ ] An ancillary program to estimate the theoretical sensitivity? Based on the extrapolate values from the alma online calculator?
- [ ] A program automatically find the good images. (Based on imstat results, and the fluxval if possible)
- [ ] modify the gen_image program to make images with different resolution, maybe three for each image? with resolution decrease three time one time

Another plan is making all your daily records online for the collaborators.

---

**2020-12-10**

Todo

- [ ] Analysis the most prominent ear in J0217-0820
- [ ] check the direction of ear with the direction of psf

Possible ways to improve the point source subtraction:

- selfcal with `nterms=2`   
- selfcal with multiscale algorithms

-----------