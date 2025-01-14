////* Code adapted: https://www.w3schools.com/howto/howto_css_modals.asp     *////

function ModalFunc(modal_id, button_id, close_id) {
    // Get the modal
    // When the user clicks on the button, open the modal
    var modal = document.getElementById(modal_id);
    modal.style.display = "block";

    // Get the button that opens the modal
    var btn = document.getElementById(button_id);

    // Get the <span> element that closes the modal
    var span = document.getElementsByClassName("close-modal")[close_id];

    // When the user clicks on <span> (x), close the modal
    span.onclick = function() {
        modal.style.display = "none";
    }

    // When the user clicks anywhere outside of the modal, close it
    window.onclick = function(event) {
        if (event.target == modal) {
            modal.style.display = "none";
        }
    }
}