/* shared quiz engine: multiple-choice self-tests per lesson.
   markup: <div class="quiz"> containing .q blocks with data-answer="x",
   radio inputs value a/b/c/d, a check button, and a .quiz-result div. */
function checkQuiz(btn) {
  var quiz = btn.closest('.quiz');
  var qs = quiz.querySelectorAll('.q');
  var right = 0, total = qs.length;
  qs.forEach(function (q) {
    var ans = q.getAttribute('data-answer');
    var sel = q.querySelector('input:checked');
    q.classList.remove('q-right', 'q-wrong');
    if (sel && sel.value === ans) {
      right++; q.classList.add('q-right');
    } else {
      q.classList.add('q-wrong');
    }
  });
  var res = quiz.querySelector('.quiz-result');
  var lang = document.documentElement.getAttribute('data-lang');
  var msg;
  if (right === total) {
    msg = lang === 'en' ? 'All ' + total + ' correct · well done!' : 'Alle ' + total + ' goed · uitstekend!';
  } else {
    msg = lang === 'en'
      ? right + ' of ' + total + ' correct · the highlighted questions deserve another look.'
      : right + ' van ' + total + ' goed · de gemarkeerde vragen verdienen nog een blik.';
  }
  res.textContent = msg;
  res.style.display = 'block';
}
