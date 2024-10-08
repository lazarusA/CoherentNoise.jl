import{_ as e,c as a,a5 as i,o as t}from"./chunks/framework.DLCQLSHW.js";const k=JSON.parse('{"title":"","description":"","frontmatter":{},"headers":[],"relativePath":"install.md","filePath":"install.md","lastUpdated":null}'),n={name:"install.md"};function l(o,s,p,r,h,d){return t(),a("div",null,s[0]||(s[0]=[i(`<h2 id="installation" tabindex="-1">Installation <a class="header-anchor" href="#installation" aria-label="Permalink to &quot;Installation {#installation}&quot;">​</a></h2><p>CoherentNoise.jl exists in the Julia&#39;s General registry. In Julia ≥ 1.0, you can add it using the package manager with:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> add CoherentNoise</span></span></code></pre></div><h2 id="Basic-Usage" tabindex="-1">Basic Usage <a class="header-anchor" href="#Basic-Usage" aria-label="Permalink to &quot;Basic Usage {#Basic-Usage}&quot;">​</a></h2><p>The following is only a brief explanation of basic usage. For more advanced usage examples, please see the <a href="/CoherentNoise.jl/dev/tutorial#first_steps">Tutorial</a> section.</p><p>To get started, let&#39;s get a feel for how to create a sampler for one of the supported noise algorithms, and sample from it. Perlin Improved noise is a well-known algorithm, so let&#39;s create a 2-dimensional Perlin Improved noise sampler and sample from it:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CoherentNoise</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sampler </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> perlin_2d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sampler, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">120.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">42.8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>0.27306901504000125</span></span></code></pre></div><p>To create a sampler, we simply call the function corresponding to the noise algorithm and dimensionality. We can then sample from it using the <code>sample</code> function, passing in as arguments a sampler, and multiple <code>Real</code> arguments corresponding to the coordinates in the sampler&#39;s noise space. In this example, we have a 2-dimensional sampler, so we passed 2 numbers; one for the coordinate along the X axis, and one for the coordinate along the Y axis.</p><div class="tip custom-block"><p class="custom-block-title">Note</p><p>It is strongly recommended to use floating-point numbers with a fractional component (non-zero value after the decimal point), as some algorithms (notably older gradient noises like Perlin Noise or Perlin &quot;Improved Noise&quot;), will return zero for integral coordinates.</p></div><p>All samplers have their own distinct random number generator, and it is how we can retrieve different results with the same input coordinates. By default, a sampler&#39;s seed is set to <code>nothing</code>, which corresponds to using your machine&#39;s hardware random number generator to create a seed. This will result in non-deterministic results, as each instance of a sampler will produce different results when sampled from. We can change the seed on a per-sampler basis, by passing the <code>seed</code> keyword argument, in order to have reproducible results:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sampler </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> perlin_2d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(seed</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">42</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sample</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sampler, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">120.2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">42.8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>-0.09267200000000297</span></span></code></pre></div><p>In summary, one creates a sampler by calling a function for the desired algorithm and dimensionality, and we can alter the the sampler&#39;s output by changing its seed.</p><p>Of particular note is that all samplers accept a <code>seed</code> keyword argument; even those that don&#39;t make use of any random numbers in their implementation. This is required for the composition pipeline feature described in the <a href="/CoherentNoise.jl/dev/tutorial#first_steps">Tutorial</a>.</p>`,15)]))}const g=e(n,[["render",l]]);export{k as __pageData,g as default};
