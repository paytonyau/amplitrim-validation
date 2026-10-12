(function () {
  var root = document.documentElement, list = document.querySelector('header.site nav ul');
  if (!list) return;
  var button = document.createElement('button');
  button.type = 'button';
  button.className = 'theme';
  button.id = 'theme-toggle';
  var item = document.createElement('li');
  item.appendChild(button);
  list.appendChild(item);
  var query = window.matchMedia ? window.matchMedia('(prefers-color-scheme: dark)') : null;
  function current() {
    var chosen = root.getAttribute('data-theme');
    if (chosen) return chosen;
    return query && query.matches ? 'dark' : 'light';
  }
  function show() {
    var dark = current() === 'dark';
    button.textContent = dark ? '\u2600 Light' : '\u263E Dark';
    button.setAttribute('aria-label', dark ? 'Switch to the light theme' : 'Switch to the dark theme');
  }
  button.addEventListener('click', function () {
    var next = current() === 'dark' ? 'light' : 'dark';
    root.setAttribute('data-theme', next);
    try { localStorage.setItem('amplitrim-theme', next); } catch (e) {}
    show();
  });
  if (query && query.addEventListener) query.addEventListener('change', show);
  show();
})();
