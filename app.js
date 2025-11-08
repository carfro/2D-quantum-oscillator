/* Lightweight OrbitControls shim (mouse rotate + zoom) for file:// usage */
(() => {
  if (typeof THREE === 'undefined' || THREE.OrbitControls) return;
  // --- minimal OrbitControls (rotate LMB, zoom wheel/middle) ---
  const STATE = { NONE: -1, ROTATE: 0, DOLLY: 1, TOUCH_ROTATE: 2, TOUCH_DOLLY_PAN: 3 };
  function OrbitControls(object, domElement) {
    if (!domElement) throw new Error('DOM element required.');
    this.object = object; this.domElement = domElement;
    this.enableDamping = true; this.dampingFactor = 0.08;
    this.enableZoom = true; this.zoomSpeed = 0.9;
    this.enableRotate = true; this.rotateSpeed = 0.9;
    this.enablePan = false; this.target = new THREE.Vector3();
    const sph = new THREE.Spherical(), dS = new THREE.Spherical();
    const q = new THREE.Quaternion().setFromUnitVectors(object.up, new THREE.Vector3(0,1,0));
    const qi = q.clone().invert(); let state = STATE.NONE, scale = 1;
    const rs = new THREE.Vector2(), re = new THREE.Vector2(), rd = new THREE.Vector2();
    const ds = new THREE.Vector2(), de = new THREE.Vector2(), dd = new THREE.Vector2();
    const pointers = []; const pointerPositions = {};
    const getSecondPointerPosition = (e) => {
      const pointer = (e.pointerId === pointers[0]) ? pointers[1] : pointers[0];
      return pointerPositions[pointer];
    };
    const onPD = (e)=>{
      pointers.push(e.pointerId);
      pointerPositions[e.pointerId] = { x: e.clientX, y: e.clientY };
      if(pointers.length === 1){
        if(e.pointerType === 'touch'){
          rs.set(e.clientX,e.clientY); state=STATE.TOUCH_ROTATE;
        } else if(e.button===0){
          rs.set(e.clientX,e.clientY); state=STATE.ROTATE;
        } else if(e.button===1){
          ds.set(e.clientX,e.clientY); state=STATE.DOLLY;
        }
      } else if(pointers.length === 2 && e.pointerType === 'touch'){
        const p2 = getSecondPointerPosition(e);
        const dx = e.clientX - p2.x, dy = e.clientY - p2.y;
        ds.set(0, Math.sqrt(dx*dx + dy*dy)); state=STATE.TOUCH_DOLLY_PAN;
      }
      if(pointers.length <= 2) domElement.setPointerCapture(e.pointerId);
    };
    const onPM = (e)=>{
      if(pointers.length === 1 && state===STATE.ROTATE){
        re.set(e.clientX,e.clientY); rd.subVectors(re,rs).multiplyScalar(this.rotateSpeed);
        const el = domElement.clientHeight;
        dS.theta -= (2*Math.PI*rd.x)/el; dS.phi -= (2*Math.PI*rd.y)/el; rs.copy(re);
        this.update();
      } else if(pointers.length === 1 && state===STATE.TOUCH_ROTATE){
        re.set(e.clientX,e.clientY); rd.subVectors(re,rs).multiplyScalar(this.rotateSpeed);
        const el = domElement.clientHeight;
        dS.theta -= (2*Math.PI*rd.x)/el; dS.phi -= (2*Math.PI*rd.y)/el; rs.copy(re);
        this.update();
      } else if(state===STATE.DOLLY){
        de.set(e.clientX,e.clientY); dd.subVectors(de,ds);
        scale *= (dd.y>0? 1/Math.pow(0.95,this.zoomSpeed): Math.pow(0.95,this.zoomSpeed)); ds.copy(de);
        this.update();
      } else if(pointers.length === 2 && state===STATE.TOUCH_DOLLY_PAN){
        pointerPositions[e.pointerId] = { x: e.clientX, y: e.clientY };
        const p2 = getSecondPointerPosition(e);
        const dx = e.clientX - p2.x, dy = e.clientY - p2.y;
        const distance = Math.sqrt(dx*dx + dy*dy);
        de.set(0, distance); dd.subVectors(de,ds);
        scale *= (dd.y>0? Math.pow(0.95,this.zoomSpeed): 1/Math.pow(0.95,this.zoomSpeed)); ds.copy(de);
        this.update();
      }
    };
    const onPU = (e)=>{
      const idx = pointers.indexOf(e.pointerId);
      if(idx !== -1) pointers.splice(idx, 1);
      delete pointerPositions[e.pointerId];
      domElement.releasePointerCapture(e.pointerId);
      if(pointers.length === 0) state=STATE.NONE;
    };
    const onWheel=(e)=>{ e.preventDefault(); scale *= (e.deltaY>0? 1/Math.pow(0.95,this.zoomSpeed): Math.pow(0.95,this.zoomSpeed)); this.update(); };
    this.update = (()=>{ const off = new THREE.Vector3(); return ()=>{
      const pos=this.object.position; off.copy(pos).sub(this.target).applyQuaternion(q);
      sph.setFromVector3(off); sph.theta += dS.theta; sph.phi += dS.phi; sph.makeSafe();
      sph.radius *= scale; this.target.addScaledVector(new THREE.Vector3(),0.0);
      off.setFromSpherical(sph).applyQuaternion(qi); pos.copy(this.target).add(off);
      this.object.lookAt(this.target); if(this.enableDamping){ dS.theta*=1-this.dampingFactor; dS.phi*=1-this.dampingFactor; } else { dS.set(0,0,0); }
      scale=1;
    };})();
    this.reset = ()=>{ this.target.set(0,0,0); this.object.position.set(0,0,11); this.object.lookAt(this.target); this.update(); };
    domElement.addEventListener('pointerdown',onPD); domElement.addEventListener('pointermove',onPM);
    domElement.addEventListener('pointerup',onPU); domElement.addEventListener('wheel',onWheel,{passive:false});
    this.update();
  }
  THREE.OrbitControls = OrbitControls;
})();

/* Quantum Harmonic Oscillator — Hermite–Gaussian 2D
 * Rendered as a dense point cloud with *spherical* sprites (lit impostors)
 * to match the exact “small balls” look.
 */
(() => {
  const wrap = document.getElementById('canvas-wrap');

  // --- Three.js scene ---
  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setSize(wrap.clientWidth, wrap.clientHeight);
  renderer.setClearColor(0x0b0c10, 1);
  wrap.appendChild(renderer.domElement);

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(42, wrap.clientWidth / wrap.clientHeight, 0.1, 100);
  camera.position.set(0.0, 0.0, 11.0);

  const DOMAIN = 4.0;       // x,y in [-DOMAIN, DOMAIN]
  const RES = 160;          // reduced grid for giant pearls (~26k points)
  const STEP = (DOMAIN * 2.0) / (RES - 1);
  const HALF_X = DOMAIN;
  const HALF_Y = DOMAIN;
  const HALF_Z = 1.0;
  const fitPadding = 1.08; // slight margin so dots do not clip screen edges

  const controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.enablePan = false;
  controls.rotateSpeed = 0.7;
  controls.zoomSpeed = 0.6;
  controls.target.set(0, 0, 0);

  const clock = new THREE.Clock();
  const defaultView = {
    position: new THREE.Vector3(0, 0, camera.position.z),
    target: new THREE.Vector3(0, 0, 0)
  };
  const cameraTween = {
    active: false,
    fromPos: new THREE.Vector3(),
    toPos: new THREE.Vector3(),
    fromTarget: new THREE.Vector3(),
    toTarget: new THREE.Vector3(),
    start: 0,
    duration: 0.8
  };

  function startCameraTween(destPos, destTarget, duration = 0.8) {
    if (camera.position.distanceTo(destPos) < 1e-6 &&
        controls.target.distanceTo(destTarget) < 1e-6) {
      cameraTween.active = false;
      return;
    }
    cameraTween.fromPos.copy(camera.position);
    cameraTween.fromTarget.copy(controls.target);
    cameraTween.toPos.copy(destPos);
    cameraTween.toTarget.copy(destTarget);
    cameraTween.start = clock.getElapsedTime();
    cameraTween.duration = duration;
    cameraTween.active = true;
  }

  function fitView(opts = {}) {
    const { animate = false } = opts;
    const halfFov = THREE.MathUtils.degToRad(camera.fov * 0.5);
    const tanHalfFov = Math.tan(halfFov);
    const aspect = renderer.domElement.clientWidth / Math.max(1, renderer.domElement.clientHeight);

    const distY = HALF_Y / tanHalfFov;
    const distX = HALF_X / (tanHalfFov * aspect);
    const dist = Math.max(distY, distX, HALF_Z * fitPadding) * fitPadding;

    camera.near = 0.1;
    camera.far = Math.max(100, dist + HALF_X * 4);
    camera.updateProjectionMatrix();

    defaultView.position.set(0, 0, dist);
    defaultView.target.set(0, 0, 0);

    if (animate) {
      startCameraTween(defaultView.position, defaultView.target);
    } else {
      cameraTween.active = false;
      camera.position.copy(defaultView.position);
      controls.target.copy(defaultView.target);
      camera.lookAt(defaultView.target);
      controls.update();
    }

    controls.target0 = defaultView.target.clone();
    controls.position0 = defaultView.position.clone();
    controls.zoom0 = camera.zoom;
  }

  // --- Geometry (regular grid packed as points) ---
  const COUNT = RES * RES;

  const geom = new THREE.BufferGeometry();
  const sampleX = new Float32Array(RES);
  const sampleY = new Float32Array(RES);
  const aXY = new Float32Array(COUNT * 2);
  let ptr = 0;
  for (let j = 0; j < RES; j++) {
    const v = j / (RES - 1);
    const y = (v * 2.0 - 1.0) * DOMAIN;
    sampleY[j] = y;
    for (let i = 0; i < RES; i++) {
      const u = i / (RES - 1);
      const x = (u * 2.0 - 1.0) * DOMAIN;
      if (j === 0) sampleX[i] = x;
      aXY[ptr++] = x;
      aXY[ptr++] = y;
    }
  }
  geom.setAttribute('position', new THREE.Float32BufferAttribute(new Float32Array(COUNT * 3), 3));
  geom.setAttribute('aXY', new THREE.BufferAttribute(aXY, 2));
  geom.boundingSphere = new THREE.Sphere(new THREE.Vector3(0, 0, 0), DOMAIN * Math.SQRT2 + 6.0);

  const SQRT2 = Math.SQRT2;
  const INV_PI_QUART = 0.7511255444649425;
  const phiCacheX = new Map();
  const phiCacheY = new Map();

  function phi1D(x, n) {
    const g = Math.exp(-0.5 * x * x);
    let phi0 = INV_PI_QUART * g;
    if (n === 0) return phi0;
    let phi1 = SQRT2 * x * phi0;
    if (n === 1) return phi1;
    let pm1 = phi0;
    let p = phi1;
    for (let i = 1; i < n; ++i) {
      const next = (SQRT2 * x * p - Math.sqrt(i) * pm1) / Math.sqrt(i + 1);
      pm1 = p;
      p = next;
    }
    return p;
  }

  function getPhi(cache, samples, n) {
    if (!cache.has(n)) {
      const arr = new Float32Array(samples.length);
      for (let i = 0; i < samples.length; ++i) {
        arr[i] = phi1D(samples[i], n);
      }
      cache.set(n, arr);
    }
    return cache.get(n);
  }

  function computeZStats(n0, n1, mix, amp) {
    const phiX0 = getPhi(phiCacheX, sampleX, n0);
    const phiX1 = getPhi(phiCacheX, sampleX, n1);
    const phiY0 = getPhi(phiCacheY, sampleY, n0);
    const phiY1 = getPhi(phiCacheY, sampleY, n1);

    let min = Infinity;
    let max = -Infinity;
    let maxAbs = 0;
    for (let j = 0; j < RES; ++j) {
      const y0 = phiY0[j];
      const y1 = phiY1[j];
      for (let i = 0; i < RES; ++i) {
        const psi0 = phiX0[i] * y0;
        const psi1 = phiX1[i] * y1;
        const psi = psi0 + (psi1 - psi0) * mix;
        const z = psi * amp;
        if (z < min) min = z;
        if (z > max) max = z;
        const absZ = Math.abs(z);
        if (absZ > maxAbs) maxAbs = absZ;
      }
    }
    return { min, max, maxAbs };
  }

  // --- Shaders ---
  const vertexShader = `
    precision highp float;
    precision highp int;
    #define MAX_N 64

    attribute vec2 aXY;

    uniform int u_n0;
    uniform int u_n1;
    uniform float u_mix;
    uniform float u_amp;
    uniform float u_pointSize;
    uniform float u_dpr;
    uniform float u_dotPx;
    uniform float u_jitter;
    uniform float u_zScale;

    varying float vI;
    varying float vZ;
    varying float vDist;

    // hash utility for blue-noise style jitter
    float hash21(vec2 p) {
      return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453123);
    }

    // normalized 1D Hermite–Gaussian via stable three-term recurrence
    float phi_n(float x, int n) {
      float g = exp(-0.5 * x * x);
      float phi0 = 0.7511255444649425 * g;           // π^{-1/4}
      if (n == 0) return phi0;
      float phi1 = 1.4142135623730951 * x * phi0;    // √2 x φ0
      if (n == 1) return phi1;
      float pm1 = phi0;
      float p   = phi1;
      for (int i = 1; i < MAX_N; ++i) {
        if (i >= n) break;
        float fi = float(i);
        float next = (1.4142135623730951 * x * p - sqrt(fi) * pm1) / sqrt(fi + 1.0);
        pm1 = p; p = next;
      }
      return p;
    }

    void main() {
      float x = aXY.x;
      float y = aXY.y;

      float phx0 = phi_n(x, u_n0);
      float phy0 = phi_n(y, u_n0);
      float psi0 = phx0 * phy0;

      float phx1 = phi_n(x, u_n1);
      float phy1 = phi_n(y, u_n1);
      float psi1 = phx1 * phy1;

      float psi = mix(psi0, psi1, u_mix);
      float zRaw = psi * u_amp;
      float zNorm = u_zScale > 0.0 ? zRaw * u_zScale : 0.0;

      vec2 jitter = (vec2(hash21(vec2(x, y)), hash21(vec2(y, x))) - 0.5) * u_jitter;
      vec3 pos = vec3(x + jitter.x, y + jitter.y, zNorm);

      vec4 mv = modelViewMatrix * vec4(pos, 1.0);
      gl_Position = projectionMatrix * mv;

      float dist = max(0.0001, length(mv.xyz));
      float px = u_dotPx * u_dpr;
      gl_PointSize = clamp(px * u_pointSize / dist, 1.0, 10.0 * u_dpr);

      // intensity used for brightness/alpha
      vI = abs(psi);
      vZ = clamp(zNorm, -1.0, 1.0);
      vDist = dist;
    }
  `;

  // NOTE: fragment shader draws a *lit sphere* in each point using gl_PointCoord
  const fragmentShader = `
    precision highp float;
    varying float vI;
    varying float vZ;
    varying float vDist;

    float srgbChannelToLinear(float c) {
      return (c <= 0.04045) ? (c / 12.92) : pow((c + 0.055) / 1.055, 2.4);
    }
    vec3 srgbToLinear(vec3 c) {
      return vec3(
        srgbChannelToLinear(c.r),
        srgbChannelToLinear(c.g),
        srgbChannelToLinear(c.b)
      );
    }
    float linearChannelToSrgb(float c) {
      return (c <= 0.0031308) ? (c * 12.92) : 1.055 * pow(c, 1.0 / 2.4) - 0.055;
    }
    vec3 linearToSrgb(vec3 c) {
      return vec3(
        linearChannelToSrgb(c.r),
        linearChannelToSrgb(c.g),
        linearChannelToSrgb(c.b)
      );
    }

    void main() {
      // map to [-1,1]^2
      vec2 p = gl_PointCoord * 2.0 - 1.0;
      float r2 = dot(p, p);
      if (r2 > 1.0) discard;          // circular sprite

      // reconstruct sphere normal from screen-space coords
      float z = sqrt(1.0 - r2);
      vec3 n = normalize(vec3(p, z));

      // single view-space light for “pearl” look
      vec3 L = normalize(vec3(0.35, 0.55, 0.76));
      float lambert = max(dot(n, L), 0.0);

      // subtle rim so dots pop on dark bg
      float rim = pow(1.0 - max(n.z, 0.0), 1.6) * 0.25;

      // contrast curve from ψ amplitude
      float gamma = 0.65;
      float I = pow(vI, gamma);
      float shadeBase = 0.42 + 0.58 * lambert + rim;

      // height-based remaps:
      // crest     – normalized height in [0,1]
      // fadeWide  – extended falloff controlling blue→white ramp
      // crestFocus/crestGlow – accentuate the very top of the surface
      float crest = clamp((vZ + 1.0) * 0.5, 0.0, 1.0);
      float fadeWide = smoothstep(-0.9, 1.2, vZ);
      float crestFocus = smoothstep(0.18, 0.95, crest);
      float crestGlow = smoothstep(0.48, 0.95, crest);

      // Diverging color map in linear RGB: negative trough (dark blue) → zero (neutral grey) → positive crest (bright white)
      vec3 negColorLin = srgbToLinear(vec3(0.48, 0.58, 0.94));
      vec3 zeroColorLin = srgbToLinear(vec3(0.40, 0.42, 0.46));
      vec3 posColorLin = srgbToLinear(vec3(1.03, 1.02, 0.98));

      float posWeight = smoothstep(0.45, 1.0, vZ);
      float negWeight = smoothstep(-1.0, -0.45, vZ);
      // Blend toward the appropriate endpoint, then return to display (sRGB) space
      vec3 baseLin = mix(negColorLin, zeroColorLin, negWeight);
      baseLin = mix(baseLin, posColorLin, posWeight);
      vec3 baseTone = linearToSrgb(baseLin);

      // gentle warm-to-cool tint
      vec3 tintLow = vec3(0.56, 0.66, 0.90);
      vec3 tintHigh = vec3(1.08, 1.04, 0.94);
      vec3 tint = mix(tintLow, tintHigh, crestFocus);

      float shadeStretch = mix(0.55, 1.38, fadeWide);
      float crestBoost = mix(1.0, 3.6, crestFocus);
      float apexGlow = mix(0.0, 2.6, crestGlow);
      float densityGain = clamp(1.0 + 2.2 / (vDist + 0.35), 1.0, 4.8);

      vec3 col = baseTone * tint * shadeBase * shadeStretch * (0.55 + 0.22 * I);
      col *= crestBoost;

      vec3 glow = tintHigh * apexGlow;
      vec3 finalCol = (col + glow) * densityGain;

      // float alpha = mix(0.36, 0.84, fadeWide) * (0.56 + 0.24 * I) * smoothstep(1.0, 0.82, r2);
      float alpha = 0.64 * smoothstep(1.0, 0.82, r2);
      alpha = min(alpha * densityGain, 1.0);

      gl_FragColor = vec4(finalCol, alpha);
      // additive blend set on material => glow builds up like the reference
    }
  `;

  const material = new THREE.ShaderMaterial({
    vertexShader,
    fragmentShader,
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
    dithering: true,
    uniforms: {
      u_n0:        { value: 0 },
      u_n1:        { value: 0 },
      u_mix:       { value: 1 },
      u_amp:       { value: 2.6 },     // vertical scale
      u_pointSize: { value: 1.0 },     // breathing multiplier
      u_dpr:       { value: Math.min(window.devicePixelRatio || 1, 2) },
      u_dotPx:     { value: 8.0 },     // target dot radius in px (doubled)
      u_zScale:    { value: 1 },
      u_jitter:    { value: STEP * 0.35 }
    }
  });

  const points = new THREE.Points(geom, material);
  scene.add(points);
  let currentTarget = 0;
  const animation = { active: false, from: 0, to: 0, start: 0, duration: 0.65 };

  function completeAnimation() {
    if (!animation.active) return;
    material.uniforms.u_n0.value = animation.to;
    material.uniforms.u_n1.value = animation.to;
    material.uniforms.u_mix.value = 1;
    animation.active = false;
    animation.from = animation.to;
    currentTarget = animation.to;
  }

  function startMorph(target) {
    if (target === currentTarget && !animation.active) return;
    if (animation.active) completeAnimation();
    if (target === currentTarget) return;

    material.uniforms.u_n0.value = currentTarget;
    material.uniforms.u_n1.value = target;
    material.uniforms.u_mix.value = 0;

    animation.from = currentTarget;
    animation.to = target;
    animation.start = clock.getElapsedTime();
    animation.active = true;
    currentTarget = target;
  }

  // --- UI wiring ---
  const quantumHeading = document.getElementById('quantumHeading');
  const nSlider = document.getElementById('nSlider');
  const sliderRow = nSlider.closest('.controls__sliderRow');
  const sliderTicks = sliderRow ? sliderRow.querySelector('.slider__ticks') : null;
  const stateEl = document.getElementById('statePsi');
  const energyE = document.getElementById('energyE');
  const degenEl = document.getElementById('degeneracyG');
  const resetBtn= document.getElementById('resetBtn');
  const theoryPanel = document.getElementById('theoryPanel');
  const theorySummary = theoryPanel ? theoryPanel.querySelector('.theory__summary') : null;
  const theoryBody = theoryPanel ? theoryPanel.querySelector('.theory__body') : null;
  const theorySummaryLabel = theorySummary ? theorySummary.querySelector('.theory__summaryLabel') : null;
  const theoryHeading = theoryBody ? theoryBody.querySelector('.theory__heading') : null;
  const theoryHeadingPrimary = theoryHeading ? theoryHeading.querySelector('.theory__headingPrimary') : null;
  const theoryHeadingSecondary = theoryHeading ? theoryHeading.querySelector('.theory__headingSecondary') : null;
  const sliderMin = parseFloat(nSlider.min || '0');
  const sliderMax = parseFloat(nSlider.max || '12');
  const sliderStep = parseFloat(nSlider.step || '1');
  const sliderSpan = sliderMax - sliderMin;
  const sliderSteps = Math.round(sliderSpan / sliderStep);
  const tickSpans = [];
  const trackMetrics = { left: 0, width: 0 };

  if (sliderTicks && sliderRow) {
    sliderTicks.innerHTML = '';
    sliderRow.style.setProperty('--tick-count', sliderSteps + 1);
    for (let i = 0; i <= sliderSteps; i++) {
      const tick = document.createElement('span');
      sliderTicks.appendChild(tick);
      tickSpans.push(tick);
    }
  }

  layoutTicks();

  function layoutTicks() {
    if (!sliderRow) return;
    const rowRect = sliderRow.getBoundingClientRect();
    const sliderRect = nSlider.getBoundingClientRect();
    const computed = window.getComputedStyle(sliderRow);
    const thumbSize = parseFloat(computed.getPropertyValue('--slider-thumb-size')) || 20;
    const trackLeftRaw = sliderRect.left - rowRect.left;
    const trackWidthRaw = sliderRect.width;
    const trackStart = trackLeftRaw + thumbSize * 0.5;
    const effectiveWidth = Math.max(0, trackWidthRaw - thumbSize);
    const trackEnd = trackStart + effectiveWidth;

    trackMetrics.left = trackStart;
    trackMetrics.width = effectiveWidth;
    trackMetrics.thumb = thumbSize;

    sliderRow.style.setProperty('--slider-track-left', `${trackStart}px`);
    sliderRow.style.setProperty('--slider-track-right', `${Math.max(0, rowRect.width - trackEnd)}px`);
    sliderRow.style.setProperty('--slider-track-width', `${Math.max(0, effectiveWidth)}px`);
    if (sliderTicks) {
      tickSpans.forEach((tick, idx) => {
        const ratio = sliderSteps === 0 ? 0 : idx / sliderSteps;
        const x = ratio * effectiveWidth;
        tick.style.left = `${x}px`;
        // Hide edge ticks to avoid visual conflict with rounded rail caps
        tick.style.opacity = (idx === 0 || idx === sliderSteps) ? '0' : '';
      });
    }
    updateSliderDecor(parseFloat(nSlider.value || sliderMin));
  }

  function updateSliderDecor(value) {
    const progressRaw = sliderSpan === 0 ? 0 : (value - sliderMin) / sliderSpan;
    const progress = Math.min(1, Math.max(0, progressRaw));
    if (sliderRow) {
      const fillPx = Math.max(0, Math.min(trackMetrics.width, progress * trackMetrics.width));
      sliderRow.style.setProperty('--slider-progress-px', `${fillPx}px`);
    }
    if (!tickSpans.length) return;
    tickSpans.forEach((tick, idx) => {
      const ratio = sliderSteps === 0 ? 0 : idx / sliderSteps;
      tick.classList.toggle('is-active', ratio <= progress + 1e-6);
    });
  }

  function renderMath(el, latex) {
    el.innerHTML = `\\(${latex}\\)`;
    if (window.MathJax && MathJax.typesetPromise) {
      MathJax.typesetPromise([el]).catch(() => {});
    }
  }

  function renderHeading(n) {
    quantumHeading.innerHTML = `
      Quantum Numbers <span class="math">\\( n_x, n_y : n_x = n_y = ${n} \\)</span>
    `;
    if (window.MathJax && MathJax.typesetPromise) {
      MathJax.typesetPromise([quantumHeading]).catch(() => {});
    }
  }

  const setN = (n) => {
    startMorph(n);
    renderHeading(n);
    renderMath(stateEl, `\\psi_{${n},${n}}(x, y)`);
    const coeff = (2 * n) + 1;     // E = (2n+1) ħω ; g = 2n+1
    const energyStr = coeff === 1 ? `E = \\hbar\\omega` : `E = ${coeff}\\,\\hbar\\omega`;
    renderMath(energyE, energyStr);
    renderMath(degenEl, `${coeff}`);
    updateSliderDecor(n);
  };

  nSlider.addEventListener('input', (e) => {
    const n = parseInt(e.target.value, 10);
    setN(n);
  });

  resetBtn.addEventListener('click', () => {
    nSlider.value = 0;
    setN(0);
    fitView({ animate: true });
  });

  if (theoryPanel && theorySummary && theorySummaryLabel && theoryBody && theoryHeadingPrimary) {
    let theoryAnimating = false;

    const easeOut = 'cubic-bezier(0.2, 0.0, 0, 1)';
    const panelEaseOpen = 'cubic-bezier(0.18, 0.0, 0.0, 1)';
    const panelEaseClose = 'cubic-bezier(0.33, 0.0, 0.07, 1)';
    const panelDurations = { open: 620, close: 540 };

    const parsePx = (value) => (value ? parseFloat(value) : 0) || 0;

    const waitForAnimation = (anim, cleanup) => {
      if (!anim) {
        if (typeof cleanup === 'function') cleanup();
        return Promise.resolve();
      }
      return new Promise((resolve) => {
        const finish = () => {
          if (typeof cleanup === 'function') cleanup();
          anim.removeEventListener('finish', finish);
          anim.removeEventListener('cancel', finish);
          resolve();
        };
        anim.addEventListener('finish', finish, { once: true });
        anim.addEventListener('cancel', finish, { once: true });
      });
    };

    const playPanelFlip = ({ direction, fromRect, toRect }) => {
      const translateX = toRect.left - fromRect.left;
      const translateY = toRect.top - fromRect.top;
      const scaleX = fromRect.width ? toRect.width / fromRect.width : 1;
      const scaleY = fromRect.height ? toRect.height / fromRect.height : 1;
      const invScaleX = scaleX ? 1 / scaleX : 1;
      const invScaleY = scaleY ? 1 / scaleY : 1;

      const keyframes = direction === 'open'
        ? [
            { transform: `translate(${-translateX}px, ${-translateY}px) scale(${invScaleX}, ${invScaleY})` },
            { transform: 'translate(0px, 0px) scale(1, 1)' }
          ]
        : [
            { transform: 'translate(0px, 0px) scale(1, 1)' },
            { transform: `translate(${translateX}px, ${translateY}px) scale(${scaleX}, ${scaleY})` }
          ];

      theoryPanel.style.willChange = 'transform';
      const anim = theoryPanel.animate(keyframes, {
        duration: direction === 'open' ? panelDurations.open : panelDurations.close,
        easing: direction === 'open' ? panelEaseOpen : panelEaseClose,
        fill: 'none'
      });
      anim.addEventListener('finish', () => {
        theoryPanel.style.willChange = '';
      }, { once: true });
      anim.addEventListener('cancel', () => {
        theoryPanel.style.willChange = '';
      }, { once: true });
      return anim;
    };

    const synthesizeCollapsedRect = (expandedRect) => {
      const summaryRect = theorySummary.getBoundingClientRect();
      const panelStyles = window.getComputedStyle(theoryPanel);
      const borderTop = parsePx(panelStyles.borderTopWidth);
      const borderBottom = parsePx(panelStyles.borderBottomWidth);
      const collapsedHeight = summaryRect.height + borderTop + borderBottom;
      return {
        top: expandedRect.top,
        left: expandedRect.left,
        width: expandedRect.width,
        height: collapsedHeight
      };
    };

    const createMorphNode = (text, referenceEl) => {
      const node = document.createElement('span');
      node.className = 'theory__morphNode';
      const refStyle = window.getComputedStyle(referenceEl);
      node.textContent = text;
      node.style.fontFamily = refStyle.fontFamily;
      node.style.fontWeight = refStyle.fontWeight;
      node.style.fontSize = refStyle.fontSize;
      node.style.lineHeight = refStyle.lineHeight;
      node.style.letterSpacing = refStyle.letterSpacing;
      node.style.textTransform = refStyle.textTransform;
      node.style.whiteSpace = refStyle.whiteSpace;
      node.style.color = refStyle.color;
      node.style.transformOrigin = 'left top';
      node.style.opacity = '1';
      return node;
    };

    const fadeSummary = (target, options = {}) => {
      if (!theorySummaryLabel) return Promise.resolve();
      const end = Math.max(0, Math.min(1, target));
      const computed = parseFloat(window.getComputedStyle(theorySummaryLabel).opacity);
      const start = Number.isFinite(computed) ? computed : (end === 0 ? 1 : 0);
      if (Math.abs(start - end) < 0.005) {
        theorySummaryLabel.style.opacity = `${end}`;
        return Promise.resolve();
      }
      theorySummaryLabel.style.willChange = 'opacity';
      const anim = theorySummaryLabel.animate(
        [
          { opacity: start },
          { opacity: end }
        ],
        {
          duration: options.duration ?? 260,
          delay: options.delay ?? 0,
          easing: options.easing ?? easeOut,
          fill: 'forwards'
        }
      );
      return waitForAnimation(anim, () => {
        theorySummaryLabel.style.opacity = `${end}`;
        theorySummaryLabel.style.willChange = '';
      });
    };

    const runMorph = ({ fromRect, toRect, direction }) => {
      const morph = createMorphNode(theorySummaryLabel.textContent || '', theoryHeadingPrimary);
      const ratioX = toRect.width === 0 ? 1 : fromRect.width / toRect.width;
      const ratioY = toRect.height === 0 ? 1 : fromRect.height / toRect.height;

      const startScaleX = direction === 'open' ? ratioX : 1;
      const startScaleY = direction === 'open' ? ratioY : 1;
      const endScaleX = direction === 'open' ? 1 : (ratioX === 0 ? 1 : 1 / ratioX);
      const endScaleY = direction === 'open' ? 1 : (ratioY === 0 ? 1 : 1 / ratioY);

      const startTransform = `translate(${fromRect.left}px, ${fromRect.top}px) scale(${startScaleX}, ${startScaleY})`;
      const endTransform = `translate(${toRect.left}px, ${toRect.top}px) scale(${endScaleX}, ${endScaleY})`;

      morph.style.transform = startTransform;
      document.body.appendChild(morph);

      const animation = morph.animate(
        [
          { transform: startTransform, opacity: 1 },
          { transform: endTransform, opacity: 1 }
        ],
        {
          duration: direction === 'open' ? 620 : 540,
          easing: direction === 'open'
            ? 'cubic-bezier(0.18, 0.0, 0.0, 1)'
            : 'cubic-bezier(0.33, 0.0, 0.07, 1)',
          fill: 'forwards'
        }
      );

      return { animation, morph };
    };

    const animateTheoryOpen = () => {
      if (theoryAnimating || (theoryPanel.open && theoryPanel.dataset.state === 'open')) return;
      theoryAnimating = true;

      const collapsedRect = theoryPanel.getBoundingClientRect();
      theoryPanel.dataset.state = 'opening';
      theoryPanel.open = true;
      theorySummary.setAttribute('aria-expanded', 'true');

      requestAnimationFrame(() => {
        const expandedRect = theoryPanel.getBoundingClientRect();
        const headingRect = theoryHeadingPrimary.getBoundingClientRect();

        if (theoryBody) {
          theoryBody.style.opacity = '0';
          theoryBody.style.transform = 'translateY(14px)';
          theoryBody.style.pointerEvents = 'none';
        }
        if (theoryHeadingPrimary) {
          theoryHeadingPrimary.style.opacity = '0';
        }
        if (theoryHeadingSecondary) {
          theoryHeadingSecondary.style.opacity = '0';
        }

        const panelAnim = playPanelFlip({
          direction: 'open',
          fromRect: collapsedRect,
          toRect: expandedRect
        });

        const { animation: morphAnim, morph } = runMorph({
          fromRect: theorySummaryLabel.getBoundingClientRect(),
          toRect: headingRect,
          direction: 'open'
        });

        const summaryFade = fadeSummary(0, { duration: 280, easing: panelEaseOpen });
        theorySummaryLabel.setAttribute('aria-hidden', 'true');

        const morphDone = waitForAnimation(morphAnim, () => {
          morph.remove();
          if (theoryHeadingPrimary) {
            theoryHeadingPrimary.style.opacity = '1';
            requestAnimationFrame(() => {
              theoryHeadingPrimary.style.opacity = '';
            });
          }
          if (theoryBody) {
            waitForAnimation(
              theoryBody.animate(
                [
                  { opacity: 0, transform: 'translateY(14px)' },
                  { opacity: 1, transform: 'translateY(0)' }
                ],
                { duration: 320, easing: easeOut, fill: 'forwards' }
              ),
              () => {
                theoryBody.style.opacity = '';
                theoryBody.style.transform = '';
                theoryBody.style.pointerEvents = '';
              }
            );
          }
          if (theoryHeadingSecondary) {
            waitForAnimation(
              theoryHeadingSecondary.animate(
                [
                  { opacity: 0, transform: 'translateY(8px)' },
                  { opacity: 1, transform: 'translateY(0)' }
                ],
                { duration: 320, easing: easeOut, delay: 160, fill: 'forwards' }
              ),
              () => {
                theoryHeadingSecondary.style.opacity = '';
                theoryHeadingSecondary.style.transform = '';
              }
            );
          }
        });

        const panelDone = waitForAnimation(panelAnim);

        Promise.all([panelDone, morphDone, summaryFade]).then(() => {
          theoryPanel.dataset.state = 'open';
          theoryAnimating = false;
        });
      });
    };

    const animateTheoryClose = () => {
      if (theoryAnimating || !theoryPanel.open) return;
      theoryAnimating = true;

      const expandedRect = theoryPanel.getBoundingClientRect();
      const collapsedRect = synthesizeCollapsedRect(expandedRect);
      const headingRect = theoryHeadingPrimary.getBoundingClientRect();
      const summaryRect = theorySummaryLabel.getBoundingClientRect();

      theoryPanel.dataset.state = 'closing';

      if (theoryHeadingSecondary) {
        theoryHeadingSecondary.style.opacity = '0';
      }
      if (theoryHeadingPrimary) {
        theoryHeadingPrimary.style.opacity = '0';
      }
      let bodyPromise = Promise.resolve();
      if (theoryBody) {
        theoryBody.style.pointerEvents = 'none';
        bodyPromise = waitForAnimation(
          theoryBody.animate(
            [
              { opacity: 1, transform: 'translateY(0)' },
              { opacity: 0, transform: 'translateY(12px)' }
            ],
            { duration: 260, easing: easeOut, fill: 'forwards' }
          ),
          () => {
            theoryBody.style.opacity = '';
            theoryBody.style.transform = '';
            theoryBody.style.pointerEvents = '';
          }
        );
      }

      const panelAnim = playPanelFlip({
        direction: 'close',
        fromRect: expandedRect,
        toRect: collapsedRect
      });

      const { animation: morphAnim, morph } = runMorph({
        fromRect: headingRect,
        toRect: summaryRect,
        direction: 'close'
      });

      theorySummaryLabel.setAttribute('aria-hidden', 'false');
      const summaryFade = fadeSummary(1, { duration: 260, easing: panelEaseClose, delay: 80 });
      const morphDone = waitForAnimation(morphAnim, () => {
        morph.remove();
      });
      const panelDone = waitForAnimation(panelAnim);

      Promise.all([panelDone, morphDone, bodyPromise, summaryFade]).then(() => {
        theoryPanel.open = false;
        theoryPanel.dataset.state = 'closed';
        theorySummary.setAttribute('aria-expanded', 'false');
        if (theoryHeadingPrimary) {
          theoryHeadingPrimary.style.opacity = '';
        }
        if (theoryHeadingSecondary) {
          theoryHeadingSecondary.style.opacity = '';
          theoryHeadingSecondary.style.transform = '';
        }
        theoryAnimating = false;
        if (typeof theorySummary.focus === 'function') {
          theorySummary.focus({ preventScroll: true });
        }
      });
    };

    const toggleTheoryPanel = (explicit) => {
      const target = typeof explicit === 'boolean' ? explicit : !theoryPanel.open;
      if (target) {
        animateTheoryOpen();
      } else {
        animateTheoryClose();
      }
    };

    theorySummary.addEventListener('click', (event) => {
      event.preventDefault();
      toggleTheoryPanel();
    });

    theorySummary.addEventListener('keydown', (event) => {
      if (event.key === ' ' || event.key === 'Enter') {
        event.preventDefault();
        toggleTheoryPanel();
      }
    });

    theoryBody.addEventListener('click', (event) => {
      const selection = window.getSelection && window.getSelection();
      if (selection && selection.toString().length) return;
      event.preventDefault();
      toggleTheoryPanel(false);
    });

    theoryPanel.dataset.state = theoryPanel.open ? 'open' : 'closed';
  }

  setN(parseInt(nSlider.value, 10));
  updateSliderDecor(parseInt(nSlider.value, 10));
  fitView();

  // --- Resize + render loop ---
  function onResize() {
    renderer.setSize(wrap.clientWidth, wrap.clientHeight);
    camera.aspect = wrap.clientWidth / wrap.clientHeight;
    layoutTicks();
    fitView();
  }
  window.addEventListener('resize', onResize);

  (function tick() {
    const t = clock.getElapsedTime();
    if (animation.active) {
      const elapsed = t - animation.start;
      const norm = Math.min(elapsed / animation.duration, 1);
      const eased = 1.0 - Math.pow(1.0 - norm, 3.0);
      material.uniforms.u_mix.value = eased;
      if (norm >= 1) {
        material.uniforms.u_n0.value = animation.to;
        material.uniforms.u_n1.value = animation.to;
        material.uniforms.u_mix.value = 1;
        animation.active = false;
        animation.from = animation.to;
        currentTarget = animation.to;
      }
    }
    if (cameraTween.active) {
      const elapsed = t - cameraTween.start;
      const norm = Math.min(elapsed / cameraTween.duration, 1);
      const eased = 1.0 - Math.pow(1.0 - norm, 3.0);
      camera.position.lerpVectors(cameraTween.fromPos, cameraTween.toPos, eased);
      controls.target.lerpVectors(cameraTween.fromTarget, cameraTween.toTarget, eased);
      if (norm >= 1) {
        camera.position.copy(cameraTween.toPos);
        controls.target.copy(cameraTween.toTarget);
        camera.lookAt(cameraTween.toTarget);
        cameraTween.active = false;
      }
    }
    const mix = material.uniforms.u_mix.value;
    const n0 = material.uniforms.u_n0.value;
    const n1 = material.uniforms.u_n1.value;
    const amp = material.uniforms.u_amp.value;
    const range = computeZStats(n0, n1, mix, amp);
    const scale = range.maxAbs > 1e-6 ? (1.0 / range.maxAbs) : 0.0;
    material.uniforms.u_zScale.value = scale;
    // barely-perceptible breathing
    material.uniforms.u_pointSize.value = 1.0 + 0.03 * Math.sin(t * 0.8);
    controls.update();
    renderer.render(scene, camera);
    requestAnimationFrame(tick);
  })();
})();
