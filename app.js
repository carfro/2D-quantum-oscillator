/* Lightweight OrbitControls shim (mouse rotate + zoom) for file:// usage */
(() => {
  if (typeof THREE === 'undefined' || THREE.OrbitControls) return;
  // --- minimal OrbitControls (rotate LMB, zoom wheel/middle) ---
  const STATE = { NONE: -1, ROTATE: 0, DOLLY: 1 };
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
    const onPD = (e)=>{ if(e.button===0){rs.set(e.clientX,e.clientY);state=STATE.ROTATE;}
                        else if(e.button===1){ds.set(e.clientX,e.clientY);state=STATE.DOLLY;}
                        domElement.setPointerCapture(e.pointerId); };
    const onPM = (e)=>{ if(state===STATE.ROTATE){
                          re.set(e.clientX,e.clientY); rd.subVectors(re,rs).multiplyScalar(this.rotateSpeed);
                          const el = domElement.clientHeight;
                          dS.theta -= (2*Math.PI*rd.x)/el; dS.phi -= (2*Math.PI*rd.y)/el; rs.copy(re);
                        } else if(state===STATE.DOLLY){
                          de.set(e.clientX,e.clientY); dd.subVectors(de,ds);
                          scale *= (dd.y>0? 1/Math.pow(0.95,this.zoomSpeed): Math.pow(0.95,this.zoomSpeed)); ds.copy(de);
                        } this.update(); };
    const onPU = (e)=>{ domElement.releasePointerCapture(e.pointerId); state=STATE.NONE; };
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

  function fitView() {
    const halfFov = THREE.MathUtils.degToRad(camera.fov * 0.5);
    const tanHalfFov = Math.tan(halfFov);
    const aspect = renderer.domElement.clientWidth / Math.max(1, renderer.domElement.clientHeight);

    const distY = HALF_Y / tanHalfFov;
    const distX = HALF_X / (tanHalfFov * aspect);
    const dist = Math.max(distY, distX, HALF_Z * fitPadding) * fitPadding;

    camera.position.set(0, 0, dist);
    camera.near = 0.1;
    camera.far = Math.max(100, dist + HALF_X * 4);
    camera.updateProjectionMatrix();

    controls.target.set(0, 0, 0);
    if (controls.target0) {
      controls.target0.copy(controls.target);
    } else {
      controls.target0 = controls.target.clone();
    }
    if (controls.position0) {
      controls.position0.copy(camera.position);
    } else {
      controls.position0 = camera.position.clone();
    }
    controls.zoom0 = camera.zoom;
    controls.update();
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
    uniform float u_zScale;

    varying float vI;
    varying float vZ;

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

      vec3 pos = vec3(x, y, zNorm);

      vec4 mv = modelViewMatrix * vec4(pos, 1.0);
      gl_Position = projectionMatrix * mv;

      float dist = max(0.0001, length(mv.xyz));
      float px = u_dotPx * u_dpr;
      gl_PointSize = clamp(px * u_pointSize / dist, 1.0, 6.0 * u_dpr);

      // intensity used for brightness/alpha
      vI = abs(psi);
      vZ = clamp(zNorm, -1.0, 1.0);
    }
  `;

  // NOTE: fragment shader draws a *lit sphere* in each point using gl_PointCoord
  const fragmentShader = `
    precision highp float;
    varying float vI;
    varying float vZ;

    // fixed aesthetic color: soft light gray (matches screenshots)
    const vec3 BASE = vec3(0.90, 0.92, 0.96);

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
      float shade = 0.35 + 0.65 * lambert + rim;

      // height-based tonal mapping (-1 → base plane, +1 → crest)
      float h = clamp(0.5 * (vZ + 1.0), 0.0, 1.0);
      float emphasis = smoothstep(0.0, 1.0, pow(h, 0.75));
      float depthBoost = mix(0.18, 1.72, emphasis);
      shade *= depthBoost;

      // subtle warm-to-cool tint across height
      vec3 tintLow = vec3(0.56, 0.64, 0.88);
      vec3 tintHigh = vec3(1.10, 1.03, 0.92);
      vec3 tint = mix(tintLow, tintHigh, emphasis);

      // final color (no scene lights — baked shading)
      vec3 col = BASE * tint * shade * (1.55 + 0.35 * I);

      // alpha balances density without blowing out additive blend
      float alpha = mix(0.45, 0.70, emphasis) * (0.68 + 0.24 * I) * smoothstep(1.0, 0.82, r2);
      gl_FragColor = vec4(col, alpha);
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
      u_zScale:    { value: 1 }
    }
  });

  const points = new THREE.Points(geom, material);
  scene.add(points);

  const clock = new THREE.Clock();
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
  const nSlider = document.getElementById('nSlider');
  const nValue  = document.getElementById('nValue');
  const stateEl = document.getElementById('statePsi');
  const energyE = document.getElementById('energyE');
  const degenEl = document.getElementById('degeneracyG');
  const resetBtn= document.getElementById('resetBtn');
  const debugToggle = document.getElementById('debugToggle');
  const debugStats = document.getElementById('debugStats');
  let debugEnabled = false;
  let lastDebugUpdate = -Infinity;
  const lastRange = { min: 0, max: 0 };

  function renderMath(el, latex) {
    el.innerHTML = `\\(${latex}\\)`;
    if (window.MathJax && MathJax.typesetPromise) {
      MathJax.typesetPromise([el]).catch(() => {});
    }
  }

  function updateDebugStats(t, force = false) {
    if (!debugEnabled) return;
    if (!force && t - lastDebugUpdate < 0.2) return;
    lastDebugUpdate = t;
    const scale = material.uniforms.u_zScale.value;
    debugStats.textContent = `z-range raw: ${lastRange.min.toFixed(3)} → ${lastRange.max.toFixed(3)} | scale ≈ ${scale.toFixed(3)}`;
  }

  const setN = (n) => {
    startMorph(n);
    renderMath(nValue, `n_x = n_y = ${n}`);
    renderMath(stateEl, `\\psi_{${n},${n}}(x, y)`);
    const coeff = (2 * n) + 1;     // E = (2n+1) ħω ; g = 2n+1
    renderMath(energyE, `E = ${coeff}\\,\\hbar\\omega`);
    renderMath(degenEl, `${coeff}`);
  };

  nSlider.addEventListener('input', (e) => {
    const n = parseInt(e.target.value, 10);
    setN(n);
  });

  resetBtn.addEventListener('click', () => {
    nSlider.value = 0;
    setN(0);
    fitView();
  });

  debugToggle.addEventListener('change', (e) => {
    debugEnabled = e.target.checked;
    debugStats.classList.toggle('is-visible', debugEnabled);
    if (debugEnabled) {
      lastDebugUpdate = -Infinity;
      updateDebugStats(clock.getElapsedTime(), true);
    } else {
      debugStats.textContent = 'z-range: —';
    }
  });

  setN(parseInt(nSlider.value, 10));
  fitView();

  // --- Resize + render loop ---
  function onResize() {
    renderer.setSize(wrap.clientWidth, wrap.clientHeight);
    camera.aspect = wrap.clientWidth / wrap.clientHeight;
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
    const mix = material.uniforms.u_mix.value;
    const n0 = material.uniforms.u_n0.value;
    const n1 = material.uniforms.u_n1.value;
    const amp = material.uniforms.u_amp.value;
    const range = computeZStats(n0, n1, mix, amp);
    lastRange.min = range.min;
    lastRange.max = range.max;
    const scale = range.maxAbs > 1e-6 ? (1.0 / range.maxAbs) : 0.0;
    material.uniforms.u_zScale.value = scale;
    updateDebugStats(t);
    // barely-perceptible breathing
    material.uniforms.u_pointSize.value = 1.0 + 0.03 * Math.sin(t * 0.8);
    controls.update();
    renderer.render(scene, camera);
    requestAnimationFrame(tick);
  })();
})();
