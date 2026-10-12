(function () {
  var blocks = document.querySelectorAll('.copyable');
  Array.prototype.forEach.call(blocks, function (block) {
    var button = document.createElement('button');
    button.type = 'button';
    button.className = 'copy';
    button.textContent = 'Copy';
    button.setAttribute('aria-label', 'Copy these commands');
    button.addEventListener('click', function () {
      var text = block.querySelector('pre').textContent;
      var done = function () {
        button.textContent = 'Copied';
        setTimeout(function () { button.textContent = 'Copy'; }, 1500);
      };
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(done, function () {});
      } else {
        var area = document.createElement('textarea');
        area.value = text;
        document.body.appendChild(area);
        area.select();
        try { document.execCommand('copy'); done(); } catch (e) {}
        document.body.removeChild(area);
      }
    });
    block.appendChild(button);
  });
})();
