document.addEventListener("DOMContentLoaded", function () {
  const header = document.querySelector(".template-article main .page-header");
  if (!header || header.querySelector(".article-print-button")) return;

  if (!header.querySelector(".author")) {
    const author = document.createElement("h4");
    author.className = "author article-print-author";
    author.setAttribute("data-toc-skip", "");
    author.textContent = "Graeme L. Hickey";
    header.querySelector("h1").after(author);
  }

  const controls = document.createElement("div");
  controls.className = "article-print-controls";

  const button = document.createElement("button");
  button.type = "button";
  button.className = "btn btn-outline-secondary btn-sm article-print-button";
  button.textContent = "Print / Save as PDF";
  button.addEventListener("click", function () {
    window.print();
  });

  controls.append(button);
  header.append(controls);
});
