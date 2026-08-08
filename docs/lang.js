(function () {
  function setLang(l) {
    document.documentElement.setAttribute('data-lang', l);
    try { localStorage.setItem('collatz-lang', l); } catch (e) {}
    document.querySelectorAll('.langswitch button').forEach(function (b) {
      b.classList.toggle('active', b.dataset.lang === l);
    });
  }
  window.setLang = setLang;
  document.addEventListener('DOMContentLoaded', function () {
    var l = 'nl';
    try { l = localStorage.getItem('collatz-lang') || 'nl'; } catch (e) {}
    setLang(l);
    document.querySelectorAll('.langswitch button').forEach(function (b) {
      b.addEventListener('click', function () { setLang(b.dataset.lang); });
    });
  });
})();
