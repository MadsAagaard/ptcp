// ── Scroll-spy ───────────────────────────────────────────────────────────────
(function() {
  if (typeof IntersectionObserver === 'undefined') return;
  const sections = document.querySelectorAll('section');
  const links    = document.querySelectorAll('#sidebar-nav a');
  const observer = new IntersectionObserver(entries => {
    entries.forEach(e => {
      if (e.isIntersecting) {
        links.forEach(a => a.classList.remove('active'));
        const lnk = document.querySelector(`#sidebar-nav a[href="#${e.target.id}"]`);
        if (lnk) lnk.classList.add('active');
      }
    });
  }, { rootMargin: '-20% 0px -70% 0px' });
  sections.forEach(s => observer.observe(s));
})();
